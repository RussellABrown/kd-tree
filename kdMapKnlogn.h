/*
 * Copyright (c) 2015, 2021, 2023, 2025 Russell A. Brown
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * The following compilation defines are relevant.
 *
 * -D PREALLOCATE - If defined, all instances of KdNodes are allocated within a vector
 *                  instead of being allocated individually. This decreases the time
 *                  required to allocate and deallocate the KdNode instances.
 * 
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 *
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * -D MERGE_CUTOFF=n - A cutoff for switching from 1 to 2 threads to merge reference
 *                     arrays in the KdNode::mergeSort* functions (default 4096)
 * 
 * -D PARTITION_CUTOFF=n - A cutoff for switching from 1 to 2 threads to partition
 *                         reference arrays in KdTree::buildKdTree (default 4096)
 */

#ifndef KD_MAP_KNLOGN_H
#define KD_MAP_KNLOGN_H

#include "kdMapNode.h"

/* A cutoff for switching from 1 to 2 threads to partition reference arrays in KdTree::buildKdTree */
#ifndef PARTITION_CUTOFF
#define PARTITION_CUTOFF 4096
#endif

/* Forward references to all classes to support any order of compilation */
template <typename, typename>
class KdTreeDynamic;

template <typename, typename>
class KdTree;

template <typename, typename>
class KdNode;

template <typename, typename>
class MergeSort;

template <typename, typename>
class NearestNeighborHeap;

/* The KdTree class defines the k-d tree API. */
template <typename K, typename V=int>
class KdTree {
private:
  KdNode<K,V>* root = nullptr;

public:
  KdNode<K,V>* getRoot() {
    return root;
  }

public:
  bool isEmpty() {
    return (root == nullptr);
  }
  
#if defined(PREALLOCATE) && !defined(KD_MAP_DYNAMIC_H)
private:
  size_t entrySize;
  vector<uint8_t>* kdNodes = nullptr;
#endif

public:
  ~KdTree() {

    // If the KdNode instances are contained by a preallocated vector,
    // delete it; otherwise, delete the KdNode::root so that the
    // ~KdNode destructor will recursively delete all tuple arrays.
    //
    // However, do not delete the KdNode::root if KD_MAP_DYNAMIC_H
    // is defined, because the KdTreeDynamic::rebuildSubTree function
    // requires that the k-d tree persist. Instead, the ~KdTreeDynamic
    // destructor will delete KdNode::root.
    //
    // Also, if the KdNode instances are contained by a preallocated
    // vector, deletion of that vector will not delete the values set
    // of each KdNode instance, so delete those values sets via the
    // deleteValues function prior to deletion of the vector.
#ifndef KD_MAP_DYNAMIC_H
#ifdef PREALLOCATE
    deleteValues(root);
    delete kdNodes;
#else
    delete root;
#endif
#endif

  }

  /*
   * Delete the values set of each KdNode instance in the tree.
   *
   * Calling parameter:
   *
   * node - the root node of the tree
   */
private:
  void deleteValues(KdNode<K,V>* const node) {
    delete node->values;
    if (node->ltChild != nullptr) {
      deleteValues(node->ltChild);
    }
    if (node->gtChild != nullptr) {
      deleteValues(node->gtChild);
    }
  }

  /*
   * The buildKdTree function builds a k-d tree by recursively partitioning
   * the reference arrays and adding KdNodes to the tree.  These arrays
   * are permuted cyclically for successive levels of the tree in
   * order that sorting occur in the order x, y, z, w...
   *
   * Calling parameters:
   *
   * reference - a KdNode*** array to recursively sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode*** temporary array from which to copy sorted results;
   * start - start element of the reference array
   * end - end element of the reference array
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K,V>* buildKdTree(KdNode<K,V>*** const references,
                                  vector< vector<signed_size_t> > const& permutation,
                                  signed_size_t const start,
                                  signed_size_t const end,
                                  signed_size_t const maximumSubmitDepth,
                                  signed_size_t const depth) {

    KdNode<K,V>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    signed_size_t p = permutation.at(depth).at(permutation.at(0).size() - 1);

    // Get the number of dimensions.
    signed_size_t dim = permutation.at(0).size() - 2;

    // Obtain the reference array that corresponds to the most significant key.
    KdNode<K,V>** reference = references[permutation.at(depth).at(dim)];

    if (end == start) {

      // Only one reference was passed to this function, so add it to the tree.
      node = reference[end];
#ifdef KD_MAP_DYNAMIC_H
      node->height = 1;
#endif

    }
    else if (end == start + 1) {

      // Two references were passed to this function in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child.
      node = reference[start];
      node->gtChild = reference[end];
#ifdef KD_MAP_DYNAMIC_H
      node->gtChild->height = 1;
      node->height = 2;
#endif

    }
    else if (end == start + 2) {

      // Three references were passed to this function in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = reference[start + 1];
      node->ltChild = reference[start];
      node->gtChild = reference[end];
#ifdef KD_MAP_DYNAMIC_H
      node->ltChild->height = node->gtChild->height = 1;
      node->height = 2;
#endif

    }
    else if (end > start + 2) {

      // Four or more references were passed to this function, so the
      // median element of the reference array is chosen as the tuple
      // about which the other reference arrays will be partitioned
      // Avoid overflow when computing the median.
      signed_size_t const median = start + ((end - start) / 2);

      // Store the median element of the reference array in a new KdNode.
      node = reference[median];

      // Build both branches with child threads at as many levels of the tree
      // as possible.  Create the child threads as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so one thread will be used.
        // Initialize startIndex=1 so that the 'for' loop that partitions the
        // reference arrays will partition a number of arrays equal to dim.
        signed_size_t startIndex = 1;

        // If depth < dim-1, copy references[permut[dim]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < dim - 1) {
          startIndex = dim - depth;
          auto const dst = references[permutation.at(depth).at(0)];
          auto const tmp = references[permutation.at(depth).at(1)];
          for (int i = start; i <= end; ++i) {
            dst[i] = reference[i];
          }
          // Sort the lower half of references[permut[0]] with the current thread.
          // Ensure that the partition p cycles as x, y, z, w...
          signed_size_t p1 = (p + 1 < dim) ? p + 1 : 0;
          MergeSort<K,V>::mergeSortReferenceAscending(dst, tmp, start, median - 1, p1, dim,
                                                      maximumSubmitDepth, depth + 1);
          // Sort the upper half of references[permut[0]] with the current thread.
          MergeSort<K,V>::mergeSortReferenceAscending(dst, tmp, median + 1, end, p1, dim,
                                                      maximumSubmitDepth, depth + 1);
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        auto const tuple = node->tuple;
        for (signed_size_t i = startIndex; i < dim; ++i) {
          // Specify the source and destination reference arrays.
          auto const src = references[permutation.at(depth).at(i)];
          auto const dst = references[permutation.at(depth).at(i - 1)];

          // Fill the lower and upper halves of one reference array
          // in ascending order with the current thread.
          for (signed_size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
            auto const src_j = src[j];
            auto const compare = MergeSort<K,V>::superKeyCompare(src_j->tuple, tuple, p, dim);
            if (compare < 0) {
              dst[++lower] = src_j;
            }
            else if (compare > 0) {
              dst[++upper] = src_j;
            }
          }
        }

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(references, permutation, start, median - 1,
                                    maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, median + 1, end,
                                    maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so two threads will be used.
        // Initialize endIndex=0 so that the 'for' loop that partitions the
        // reference arrays will partition a number of arrays equal to dim.
        signed_size_t startIndex = 1;

        // If depth < dim-1, copy references[permut[dim]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < dim - 1) {
          startIndex = dim - depth;
          auto const dst = references[permutation.at(depth).at(0)];
          auto const tmp = references[permutation.at(depth).at(1)];
          // Copy and sort the lower half of references[permut[0]] with a child thread.
          // Ensure that the partition p cycles as x, y, z, w...
          signed_size_t p1 = (p + 1 < dim) ? p + 1 : 0;
          auto copyFuture =
            async(launch::async, [&] {
                                   for (int i = start; i <= median - 1; ++i) {
                                     dst[i] = reference[i];
                                   }
                                   MergeSort<K,V>::mergeSortReferenceAscending(dst, tmp, start, median - 1, p1,
                                                                               dim, maximumSubmitDepth, depth);
                                 });

          // Copy and sort the upper half of references[permut[0]] with the current thread.
          for (int i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
          }
          MergeSort<K,V>::mergeSortReferenceAscending(dst, tmp, median + 1, end, p1, dim, maximumSubmitDepth, depth);

          // Wait for the child thread to finish execution.
          try {
            copyFuture.get();
          }
          catch (exception const& e) {
            throw runtime_error("\n\ncaught exception for copy future in buildKdTree\n");
          }
        }

        // Create a copy of the node->tuple array so that the current thread
        // and the child thread do not contend for read access to this array.
        auto const tuple = node->tuple;
        auto const point = new K[dim];
        for (signed_size_t i = 0; i < dim; ++i) {
          point[i] = tuple[i];
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        for (signed_size_t i = startIndex; i < dim; ++i) {
          // Specify the source and destination reference arrays.
          auto const src = references[permutation.at(depth).at(i)];
          auto const dst = references[permutation.at(depth).at(i - 1)];

          // Two threads may be used to partition the reference arrays, analogous to
          // the use of two threads to merge the results for the merge sort algorithm.
          // Are there sufficient array elements to justify multi-threaded processing?
          if (end - start + 1 > PARTITION_CUTOFF)
           {
              // Yes, so fill one reference array in ascending order with a child thread.
            auto partitionFuture =
              async(launch::async, [&] {
                                    for (signed_size_t lower = start - 1, upper = median, j = start; j <= median; ++j) {
                                      auto const src_j = src[j];
                                      auto const compare = MergeSort<K,V>::superKeyCompare(src_j->tuple, point, p, dim);
                                      if (compare < 0) {
                                        dst[++lower] = src_j;
                                      }
                                      else if (compare > 0) {
                                        dst[++upper] = src_j;
                                      }
                                    }
                                  });

            // Simultaneously fill the same reference array in descending order with the current thread.
            for (signed_size_t lower = median, upper = end + 1, k = end; k > median; --k) {
              auto const src_k = src[k];
              auto const compare = MergeSort<K,V>::superKeyCompare(src_k->tuple, tuple, p, dim);
              if (compare < 0) {
                dst[--lower] = src_k;
              }
              else if (compare > 0) {
                dst[--upper] = src_k;
              }
            }

            // Wait for the child thread to finish execution.
            try {
              partitionFuture.get();
            }
            catch (exception const& e) {
              throw runtime_error("\n\ncaught exception for partition future in buildKdTree\n");
            }
          }
          else
          {
            // No, there are insufficient array elements to justify multi-threaded processing,
            // so fill the lower and upper halves of one reference array in ascending order
            // with the current thread.
            for (signed_size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
              auto const src_j = src[j];
              auto const compare = MergeSort<K,V>::superKeyCompare(src_j->tuple, tuple, p, dim);
              if (compare < 0) {
                dst[++lower] = src_j;
              }
              else if (compare > 0) {
                dst[++upper] = src_j;
              }
            }
          }
        }

        // Delete the point array.
        delete[] point;

        // Recursively build the < branch of the tree with a child thread.
        auto buildFuture = async(launch::async,
                                 buildKdTree,
                                 references,
                                 ref(permutation),
                                 start,
                                 median - 1,
                                 maximumSubmitDepth,
                                 depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, median + 1, end,
                                    maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for build future in buildKdTree\n");
        }
      }

#ifdef KD_MAP_DYNAMIC_H
      // Compute the height at this node as the recursion unwinds.
      node->height = KdTreeDynamic<K,V>::computeHeight(node);
#endif

    }
    else if (end < start) {

      // This is an illegal condition that should never occur, so test for it last.
      ostringstream buffer;
      buffer << "\n\nerror has occurred at depth = " << depth << " : end = " << end
             << "  <  start = " << start << " in buildKdTree\n";
      throw runtime_error(buffer.str());

    }

    // Return the pointer to the root of the k-d tree.
    return node;
  }

  /*
   * The swap function swaps two elements in a vector<signed_size_t>.
   *
   * Calling parameters:
   *
   * a - the vector
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  inline
  static void swap(vector<signed_size_t>& a,
                   signed_size_t const i,
                   signed_size_t const j) {
    
    signed_size_t const t = a[i];
    a[i] = a[j];
    a[j] = t;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   * This version of createKdTree is called from the KdTreeDynamic::rebuildSubtree function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode<K>*> whose KdNodes store the (x, y, z, w...) coordinates
   * dim - the number of dimensions (required when KD_MAP_DYNAMIC_H is defined)
   * coordinates - a vector of pairs that store the coordinates and their associated values
   * p - the leading dimension
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
public:
  static KdTree<K,V>* createKdTree(vector<KdNode<K,V>*> const& kdNodes,
                                   size_t const dim,
                                   signed_size_t const maximumSubmitDepth,
                                   signed_size_t const p) {

    // Create a KdTree instance.
    auto tree = new KdTree();

    // Allocate the references arrays including one additional array.
    size_t numDimensions = dim;
    KdNode<K,V>*** references = new KdNode<K,V>**[numDimensions + 1];
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      references[i] = new KdNode<K,V>*[kdNodes.size()];
    }

    // Don't allocate KdNodes instances for the pth references array
    // (where p is the leading dimension) instead of the first
    // references array to permit KdTreeDynamic::balanceSubtree
    // to build a sub-tree whose root node has a non-zero
    // partition coordinate.
    //
    // Copy pointers from the KdNode instances of the kdNodes vector.
    // These pointers will be re-ordered by the MergeSort member
    // functions and then copied back to the kdNode instances by the
    // KdNode::getKdNode function. For this case where KD_MAP_DYNAMIC_H
    // is defined, the tuples will be deallocated by the ~KdNode destructor
    // when a KdNode instance is deleted from the dynamic k-d tree.
    for (size_t i = 0; i < kdNodes.size(); ++i) {
      references[p][i] = kdNodes[i];
    }

    // Sort the pth references array using multiple threads, where p
    // is the leading dimension. Importantly, or compatibility with the
    // permutation vector initialized below, use the pth dimension as
    // the leading key of the super key. Also, only the pth references
    // array has been populated with KdNode instances.
    MergeSort<K,V>::mergeSortReferenceAscending(references[p], references[numDimensions],
                                                0, kdNodes.size() - 1,
                                                p, numDimensions, maximumSubmitDepth, 0);

    // For a dynamic k-d tree, it is unnecessary remove duplicate coordinates,
    // so only specify the end index.
    signed_size_t const end = kdNodes.size() - 1;

    // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() ).
    signed_size_t maxDepth = 1;
    signed_size_t size = kdNodes.size();
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call of
    // the buildKdTree function because that coordinate depends only on the depth of
    // recursion, so it may be pre-computed and stored in the permutation vector.
    //
    // first create an array the tracks which reference array contains each coordinate.
    // index 0 is for x, index 1 is for y ... index dim (i.e., numDimensions) is the
    // reference array that MergeSort::mergeSortReferenceAscending used as a temporary array.
    std::vector<signed_size_t> current(numDimensions + 1);
    for(size_t i = 0;  i <= numDimensions; ++i) {
      current[i] = i;
    }
    std::vector<signed_size_t> indices(numDimensions + 2);

    // Create a 2D 'permutation' vector from the 'indices' vector to specify permutation
    // of the reference arrays and of the partition coordinate, and a 'verifyPermutation'
    // vector to specify permutation of the partition coordinate for verifyKdTree.
    vector< vector<signed_size_t> > permutation(maxDepth, vector<signed_size_t>(numDimensions + 2));

    // Fill the permutation vector by calculating the permutation of the indices vector
    // and the the partition coordinate of the tuple at each depth in the tree.
    for (signed_size_t depth = 0; depth < maxDepth; ++depth) {
      signed_size_t p0 = (depth + p) % numDimensions;
      // The last entry of the indices vector contains the partition coordinate.
      indices[numDimensions + 1] = p0;
      // The penultimate entry of the indices vector specifies the source reference array.
      indices[numDimensions] = current[p0];
      // The first entry of the indices vector specifies the temporary reference array.
      indices[0] = current[numDimensions];
      size_t k = 1;
      // do the partitioning swaps
      for (size_t i = 1;  i < numDimensions; ++i) {
        size_t j = (i + p0) % numDimensions;  // this is the coordinate index following the primary
        indices[k] = current[j];  // write the reference index for that coordinate to the indices array
        swap(current, numDimensions, j);  // keep track of where the coordinate will have been have been swapped to
        ++k;
      }
      permutation[depth] = indices;  // write the index vector to he current level of the permutation vector
    }

    // Build the k-d tree with multiple threads if possible.
    tree->root = buildKdTree(references, permutation, 0, end, maximumSubmitDepth, 0);

    // Delete the references arrays but not the KdNodes instances that they point to
    // because those KdNodes instances will be deleted by the ~KdTree destructor.
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the KdTree instance.
    return tree;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   *
   * Calling parameters:
   *
   * coordinates - a vector of pairs that store the coordinates and their associated values
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * numberOfNodes - the number of nodes counted by KdNode::verifyKdTree - returned by reference
   * allocateTime, sortTime, removeTime, kdTime, verifyTime, deallocateTime - execution times returned by reference
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
public:
  static KdTree<K,V>* createKdTree(vector<pair<vector<K>,V>> const& coordinates,
                                   signed_size_t const maximumSubmitDepth,
                                   signed_size_t& numberOfNodes,
                                   double& allocateTime,
                                   double& sortTime,
                                   double& removeTime,
                                   double& kdTime,
                                   double& verifyTime,
                                   double& deallocateTime) {

    // Create a KdTree instance.
    auto tree = new KdTree();

    // Allocate the references arrays including one additional array.
    auto beginTime = steady_clock::now();
    size_t numDimensions = coordinates[0].first.size();
    KdNode<K,V>*** references = new KdNode<K,V>**[numDimensions + 1];
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      references[i] = new KdNode<K,V>*[coordinates.size()];
    }

#ifdef PREALLOCATE
    // Allocate all KdNodes instances as a single vector so that they
    // may be subsequently deleted as a single vector by the ~KdTree
    // destructor, which is faster than deleting them individually.
    //
    // Point each element of the first references array to a KdNode instance
    // that is an element of the kdNodes vector and initalize that instance.
    // KdNode::tuple is an array of 1 element that is extended
    // to dimensions elements by appending dimensions-1 elements
    // to the KdNode instance.
    //
    // Because KdNode::tuple contains one element of type K,
    // the alignment of KdNode at least as large as the
    // alignment of K. Round up all alignments to the next
    // multiple of kdNodeAlign.
    size_t const kdNodeAlign = alignof(KdNode<K,V>);
    size_t const kdNodeSize = ((sizeof(KdNode<K,V>) + kdNodeAlign - 1) / kdNodeAlign) * kdNodeAlign;
    size_t const setSize = ((sizeof(set<V>) + kdNodeAlign - 1) / kdNodeAlign) * kdNodeAlign;
    size_t const tupleSize = ((sizeof(K) * (numDimensions - 1)) / kdNodeAlign) * kdNodeAlign;
    tree->entrySize = kdNodeSize + setSize + tupleSize;
    // The following kdNodeAlign argument to new is likely redundant and requires c++17. See
    // https://stackoverflow.com/questions/15511909/does-the-alignas-specifier-work-with-new
    tree->kdNodes = new vector<uint8_t>(tree->entrySize * coordinates.size(), kdNodeAlign); // requires c++17
    for (size_t i = 0; i < coordinates.size(); ++i) {
      new(&(*(tree->kdNodes))[tree->entrySize * i]) KdNode<K,V>(coordinates, i);
      references[0][i] = reinterpret_cast<KdNode<K,V>*>(&(*(tree->kdNodes))[tree->entrySize * i]);
    }
#else
    // Allocate KdNode instances for the first references array. These
    // KdNode instances will be deallocated by the ~KdTree destructor.
    for (size_t i = 0; i < coordinates.size(); ++i) {
      references[0][i] = new KdNode<K,V>(coordinates, i);
    }
#endif

    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    allocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Sort the first references array using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
    // Also, only the first references array is populated with T arrays.
    beginTime = steady_clock::now();
    MergeSort<K,V>::mergeSortReferenceAscending(references[0], references[numDimensions],
                                                0, coordinates.size() - 1,
                                                0, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    sortTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Remove references to duplicate coordinates via one pass through the first reference array.
    beginTime = steady_clock::now();
    signed_size_t const end = KdNode<K,V>::removeDuplicates(references[0], 0, numDimensions, coordinates.size());
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    removeTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Start the timer to time building the k-d tree.
    beginTime = steady_clock::now();

    // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() ).
    signed_size_t maxDepth = 1;
    signed_size_t size = coordinates.size();
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute either the permutation of the reference array or
    // the partition coordinate upon each recursive call of the buildKdTree function
    // because both depend only on the depth of recursion, so they may be pre-computed.
    //
    // Because this vector is initialized with 0, 1, 2, 3, 0, 1, 2, 3, etc. (for
    // e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference
    // array above using 0 as the leading key of the super key.
    //
    // Begin by creating an 'indices' vector.
    vector<signed_size_t> indices(numDimensions + 2);
    for (size_t i = 0; i < indices.size() - 1; ++i) {
      indices[i] = i;
    }

    // Create a 2D 'permutation' vector from the 'indices' vector to specify permutation
    // of the reference arrays and of the partition coordinate.
    vector< vector<signed_size_t> > permutation(maxDepth, vector<signed_size_t>(numDimensions + 2));
    vector<signed_size_t> permutationVerify(maxDepth);

    // Fill the permutation vector by calculating the permutation of the indices vector
    // and the the partition coordinate of the tuple at each depth in the tree.
    for (size_t i = 0; i < permutation.size(); ++i) {
      // The last entry of the indices vector contains the partition coordinate.
      indices[numDimensions + 1] = permutationVerify[i] = i % numDimensions;
      // Swap the first and second to the last elements of the indices vector.
      swap(indices, 0, numDimensions);
      // Copy the indices vector to one row of the permutation vector.
      permutation[i] = indices;
      // Swap the third and second to the last elements of the indices vector.
      swap(indices, numDimensions - 1, numDimensions);
    }

    // Build the k-d tree with multiple threads if possible.
    tree->root = buildKdTree(references, permutation, 0, end, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    kdTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Verify the k-d tree and report the number of kdNodes.
    beginTime = steady_clock::now();
    numberOfNodes = tree->root->verifyKdTree(permutationVerify, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    verifyTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
  
    // Delete the references arrays but not the KdNodes instances that they point to
    // because those KdNodes instances will be deleted by the ~KdTree destructor.
    beginTime = steady_clock::now();
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    deallocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Return the pointer to the KdTree instance.
    return tree;
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list<KdNode<K,V>*> that is passed by reference and modified
   * queryLower - the query lower bound vector that is passed by reference and modified
   * queryUpper - the query upper bound vector that is passed by reference and modified
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K,V>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth) {

    root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth);
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list of KdNode pointers that is passed by reference and modified
   * queryLower - the query lower bound vector that is passed by reference and modified
   * queryUpper - the query upper bound vector that is passed by reference and modified
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K,V>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    vector<bool> const& enable) {
    
    root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth, enable);
  }

  /*
   * Walk the k-d tree and append to a list each KdNode that lies inside
   * the hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list of KdNode pointers that is passed by reference and modified
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   *
   * return a list of pointers to KdNodes that lie within the query hyper-rectangle.
   */
public:
  void bruteRegion(list<KdNode<K,V>*>& result,
                   vector<K>& queryLower,
                   vector<K>& queryUpper) {

    root->bruteRegion(result, queryLower, queryUpper);
  }

  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to attempt to find
  */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K,V>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors) {
    
    root->findNearestNeighbors(neighbors, query, numNeighbors);
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K,V>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                           vector<bool> const& enable) {
    
    root->findNearestNeighbors(neighbors, query, numNeighbors, enable);
  }
  
  /*
   * Verify the consistency between the nearest neighbors lists found
   * by k-d tree search and by brute force.
   *
   * Calling parameters:
   *
   * neighborsFast - a list of nearest neighbors found by k-d tree search
   * neighborsSlow - a list of nearest neighbors found by brute force.
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the lists. Hence, this function is not static.
   */
public:
  void verifyNearestNeighbors(forward_list< pair<double, KdNode<K,V>*> >& neighborsFast,
                              forward_list< pair<double, KdNode<K,V>*> >& neighborsSlow) const {
    
    root->verifyNearestNeighbors(neighborsFast, neighborsSlow);
  }

  /*
   * Sort a list of KdNode instances by increasing distance from the query point.
   *
   * Calling parameters:
   * 
   * kdList - the list of KdNode instances
   * query - a vector that contains the query point coordinates
   * maxNodes - the maximum number of nodes to maintain on the heap
   * 
   * returns: a sorted forward list of KdNode instances
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  forward_list<pair<double, KdNode<K,V>*>> sortByDistance(list<KdNode<K,V>*> const& kdList,
                                                          vector<K> const& query,
                                                          signed_size_t const& maxNodes) {

    return root->sortByDistance(kdList, query, maxNodes);
  }

#ifdef REVERSE_NEAREST_NEIGHBORS
  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to a nearest neighbors vector and to a reverse
   * nearest neighbors vector.
   *
   * Each element of the reverse nearest neighbors vector contains a list
   * of KdNodes to which the reference KdNode is a nearest neighbor.
   *
   * The concept of reverse nearest neighbors was first described by F. Korn and
   * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
   * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
   * Data, pages 201-212.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   */
public:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K,V>*> > >& nn,
                                   vector< forward_list< pair<double, KdNode<K,V>*> > >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t maximumSubmitDepth) {
    
    root->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions, numNeighbors, maximumSubmitDepth);
  }

  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to a nearest neighbors vector and to a reverse
   * nearest neighbors vector.
   *
   * Each element of the reverse nearest neighbors vector contains a list
   * of KdNodes to which the reference KdNode is a nearest neighbor.
   *
   * The concept of reverse nearest neighbors was first described by F. Korn and
   * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
   * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
   * Data, pages 201-212.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K,V>*> >* >& nn,
                                   vector< forward_list< pair<double, KdNode<K,V>*> >* >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t const maximumSubmitDepth,
                                   vector<bool> const& enable) {
    
    root->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions, numNeighbors, maximumSubmitDepth, enable);
   }

  /*
   * Verify the correctness of the reverse nearest neighbors vector.
   *
   * Calling parameter:
   *
   * nn - the nearest neighbors vector
   * rnn - the reverse nearest neighbors vector
   * size - the number of KdNode instances prior to calling KdMapNode::removeDuplicates
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vectors. Hence, this function is not static.
   */
public:
  void verifyReverseNeighbors(vector< forward_list< pair<double, KdNode<K,V>*> > >& nn,
                              vector< forward_list< pair<double, KdNode<K,V>*> > >& rnn,
                              size_t const size) {

    root->verifyReverseNeighbors(nn, rnn, size);
  }

  /*
   * Calculate the mean and standard deviation of the distances and list sizes in a vector.
   *
   * Calling parameter:
   *
   * vec - a vector of lists
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vector. Hence, this function is not static.
   */
public:
  void calculateMeanStd(vector< forward_list< pair<double, KdNode<K,V>*> > >& vec,
                        double& meanSize,
                        double& stdSize,
                        double& meanDist,
                        double& stdDist) const {

    root->calculateMeanStd(vec, meanSize, stdSize, meanDist, stdDist);
  }

  /*
   * Count the number of non-empty lists in a vector.
   *
   * Calling parameter:
   *
   * vec - a vector of lists
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vector. Hence, this function is not static.
   */
public:
  size_t nonEmptyLists(vector< forward_list< pair<double, KdNode<K,V>*> > >& vec) const {

    return root->nonEmptyLists(vec);
  }
#endif // REVERSE_NEAREST_NEIGHBORS

  /*
   * Find M nearest neighbors to the query vector via brute force
   * and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   */
public:
  void bruteNearestNeighbors(forward_list< pair<double, KdNode<K,V>*> >& neighbors,
                             vector<K> const& query,
                             signed_size_t const numNeighbors) {
    
    root->bruteNearestNeighbors(neighbors, query, numNeighbors);
  }
  
  /*
   * The verifyKdTree function walks the k-d tree and checks that the
   * children of a node are in the correct branch of that node.
   *
   * Calling parameters:
   *
   * permutation - the permutation vector
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns: a count of the number of kdNodes in the k-d tree
   */
public:
  signed_size_t verifyKdTree(signed_size_t const dim,
                             signed_size_t const maximumSubmitDepth) {

    if (root != nullptr) {
      return root->verifyKdTree(dim, maximumSubmitDepth, 0, 0);
    } else {
      return 0;
    }
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameters:
   *
   * tuple - the tuple as an array
   * dim - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuple(K const* tuple,
                  signed_size_t const dim) const {
    
    root->printTuple(tuple, dim);
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameter:
   *
   * tuple - the tuple as a vector
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuple(vector<K> const& tuple) const {
    
    root->printTuple(tuple);
  }

  /*
   * The printTuples function prints all tuples in a forward list of pairs.
   *
   * Calling parameters:
   *
   * regionList - a forward list of (double, KdNode*) pairs
   * maximumNumberOfNodesToPrint - the maximum number of KdNodes to print
   * numDimensions - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuples(forward_list<pair<double, KdNode<K,V>*>> const& regionList,
                   signed_size_t const maximumNumberOfNodesToPrint,
                   signed_size_t const numDimensions) const {
    
    root->printTuples(regionList, maximumNumberOfNodesToPrint, numDimensions);
  }

  /*
   * The printKdTree function prints the k-d tree "sideways" with the root at the ltChild.
   *
   * Calling parameters:
   *
   * dim - the number of dimensions
   * depth - the depth in the k-d tree
   */
public:
  void printKdTree(signed_size_t const dim,
                   signed_size_t const depth) const {

    root->printKdTree(dim, depth);
  }

  friend class KdTreeDynamic<K,V>;
}; // class KdTree

#endif // KD_MAP_KNLOGN_H