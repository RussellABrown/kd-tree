/*
 * Copyright (c) 2015, 2021, 2023, 2024, 2025 Russell A. Brown
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
 * -D PREALLOCATE - If defined, all instances of KdNodes and all tuple arrays are
 *                  allocated en masse within vectors instead of being allocated
 *                  individually. This decreases the time required to allocate and
 *                  deallocate the KdNode instances and tuple arrays.
 * 
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 *
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * 
 * -D MERGE_CUTOFF=n - A cutoff for switching from 1 to 2 threads to merge reference
 *                     arrays in the KdNode::mergeSort* functions (default 4096)
 * 
 * -D PARTITION_CUTOFF=n - A cutoff for switching from 1 to 2 threads to partition
 *                         reference arrays in KdTree::buildKdTree (default 4096)
 */

#ifndef KD_TREE_KNLOGN_H
#define KD_TREE_KNLOGN_H

#include "kdTreeNode.h"

/* A cutoff for switching from 1 to 2 threads to partition reference arrays in KdTree::buildKdTree */
#ifndef PARTITION_CUTOFF
#define PARTITION_CUTOFF 4096
#endif

/* The KdTree class defines the k-d tree API. */
template <typename K>
class KdTree {
private:
  KdNode<K>* root = nullptr;

#ifdef PREALLOCATE
  vector<KdNode<K>>* kdNodes;
  vector<K>* tuples;
#endif

public:
  ~KdTree() {

    // If the KdNode instances and the tuples contained by preallocated
    // vectors, delete them; otherwise, delete the root KdNode so that
    // the ~KdNode destructor will recursively delete all tuple arrays.
#ifdef PREALLOCATE
    delete kdNodes;
    delete tuples;
#else
    delete root;
#endif

  }

  /*
   * The buildKdTree function builds a k-d tree by recursively partitioning
   * the reference arrays and adding KdNodes to the tree.  These arrays
   * are permuted cyclically for successive levels of the tree in
   * order that sorting occur in the order x, y, z, w...
   *
   * Calling parameters:
   *
   * references - a K*** of references to each of the (x, y, z, w...) tuples
   * permutation - a vector<vector<signed_size_t>> that indicates permutation of the reference arrays
   * kdNodes - a vector<KdNode<K>>* if PREALLOCATE is defined
   *         - a vector<KdNode<K>*> if PREALLOCATE is undefined
   * start - start element of the reference array
   * end - end element of the reference array
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K>* buildKdTree(K*** const references,
                                vector<vector<signed_size_t>> const& permutation,
#ifndef PREALLOCATE
                                vector<KdNode<K>*> const& kdNodes,
#else
                                vector<KdNode<K>>* const kdNodes,
#endif
                                signed_size_t const start,
                                signed_size_t const end,
                                signed_size_t const maximumSubmitDepth,
                                signed_size_t const depth) {

    KdNode<K>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    signed_size_t const p = permutation[depth][permutation[0].size() - 1];

    // Get the number of dimensions.
    signed_size_t const dim = permutation[0].size() - 2;

    // Obtain the reference array that corresponds to the most significant key.
    auto const reference = references[permutation[depth][dim]];

    if (end == start) {

      // Only one reference was passed to this function, so add it to the tree.
      node = KdNode<K>::getKdNode(reference, kdNodes, end);

    }
    else if (end == start + 1) {

      // Two references were passed to this function in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child.
      node = KdNode<K>::getKdNode(reference, kdNodes, start);
      node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);

    }
    else if (end == start + 2) {

      // Three references were passed to this function in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = KdNode<K>::getKdNode(reference, kdNodes, start + 1);
      node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, start);
      node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);

    }
    else if (end > start + 2) {

      // Four or more references were passed to this function, so the
      // median element of the reference array is chosen as the tuple
      // about which the other reference arrays will be partitioned
      // Avoid overflow when computing the median.
      signed_size_t const median = start + ((end - start) / 2);

      // Store the median element of the reference array in a new KdNode.
      node = KdNode<K>::getKdNode(reference, kdNodes, median);

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
          auto const dst = references[permutation[depth][0]];
          auto const tmp = references[permutation[depth][1]];
          for (int i = start; i <= end; ++i) {
            dst[i] = reference[i];
          }
          // Sort the lower half of references[permut[0]] with the current thread.
          MergeSort<K>::mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1, dim,
                                                    maximumSubmitDepth, depth + 1);
          // Sort the upper half of references[permut[0]] with the current thread.
          MergeSort<K>::mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1, dim,
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
          auto const src = references[permutation[depth][i]];
          auto const dst = references[permutation[depth][i-1]];

          // Fill the lower and upper halves of one reference array
          // in ascending order with the current thread.
          for (signed_size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
            auto const src_j = src[j];
            auto const compare = MergeSort<K>::superKeyCompare(src_j, tuple, p, dim);
            if (compare < 0) {
              dst[++lower] = src_j;
            }
            else if (compare > 0) {
              dst[++upper] = src_j;
            }
          }
        }

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(references, permutation, kdNodes, start, median - 1,
                                    maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, kdNodes, median + 1, end,
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
          auto const dst = references[permutation[depth][0]];
          auto const tmp = references[permutation[depth][1]];
          // Copy and sort the lower half of references[permut[0]] with a child thread.
          auto copyFuture =
            async(launch::async, [&] {
                                   for (int i = start; i <= median - 1; ++i) {
                                     dst[i] = reference[i];
                                   }
                                   MergeSort<K>::mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1,
                                                                             dim, maximumSubmitDepth, depth);
                                 });

          // Copy and sort the upper half of references[permut[0]] with the current thread.
          for (int i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
          }
          MergeSort<K>::mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1,
                                                    dim, maximumSubmitDepth, depth);

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
        vector<K> point(dim);
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
          auto const src = references[permutation[depth][i]];
          auto const dst = references[permutation[depth][i-1]];

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
                                      auto const compare = MergeSort<K>::superKeyCompare(src_j, point.data(), p, dim);
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
              auto const compare = MergeSort<K>::superKeyCompare(src_k, tuple, p, dim);
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
              auto const compare = MergeSort<K>::superKeyCompare(src_j, tuple, p, dim);
              if (compare < 0) {
                dst[++lower] = src_j;
              }
              else if (compare > 0) {
                dst[++upper] = src_j;
              }
            }
          }
        }

        // Recursively build the < branch of the tree with a child thread.
        auto buildFuture = async(launch::async,
                                 buildKdTree,
                                 references,
                                 ref(permutation),
                                 ref(kdNodes),
                                 start,
                                 median - 1,
                                 maximumSubmitDepth,
                                 depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, kdNodes, median + 1, end,
                                    maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for build future in buildKdTree\n");
        }
      }

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
   * The swap function swaps two vector elements.
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
   *
   * Calling parameters:
   *
   * coordinates - a vector<vector<K>> that stores the (x, y, z, w...) coordinates
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * numberOfNodes - the number of nodes counted by KdNode::verifyKdTree - returned by reference
   * allocateTime, sortTime, removeTime, kdTime, verifyTime, deallocateTime - execution times returned by reference
   *
   * returns: a KdTree pointer to the root of the k-d tree
   */
public:
  static KdTree<K>* createKdTree(vector<vector<K>> const& coordinates,
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
    size_t numDimensions = coordinates[0].size();
    K*** references = new K**[numDimensions + 1];
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      references[i] = new K*[coordinates.size()];
    }

#ifndef PREALLOCATE
    // Allocate tuple arrays for the first references array.
    // Each tuple array will be deallocated by either the
    // KdNode::removeDuplicates function or the ~KdNode destructor.
    // Copy (x, y, z, w...) coordinates to each tuple array.
    for (size_t i = 0; i < coordinates.size(); ++i) {
      references[0][i] = new K[numDimensions];
      for (size_t j = 0; j < numDimensions; ++j) {
        references[0][i][j] = coordinates[i][j];
      }
    }
#else
    // Allocate a tuple vector that comprises all tuple arrays and
    // assign individual tuple arrays to the first references array.
    // This vector will be deallocated by the ~KdTree destructor.
    // Copy (x, y, z, w...) coordinates to each tuple array.
    tree->tuples = new vector<K>(numDimensions * coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      references[0][i] = &(*(tree->tuples)).data()[numDimensions * i];
      for (size_t j = 0; j < numDimensions; ++j) {
        references[0][i][j] = coordinates[i][j];
      }
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
    MergeSort<K>::mergeSortReferenceAscending(references[0], references[numDimensions],
                                              0, coordinates.size() - 1,
                                              0, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    sortTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Remove references to duplicate coordinates via one pass through the first reference array.
    beginTime = steady_clock::now();
    signed_size_t const end = KdNode<K>::removeDuplicates(references[0], 0, numDimensions, coordinates.size());
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    removeTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    beginTime = steady_clock::now();
#ifndef PREALLOCATE
    // Allocate the kdNodes vector of pointers to KdNode instances
    // but only for references that have not been removed.
    // Allocate a KdNode instance for each entry in the kdNodes vector.
    vector<KdNode<K>*> kdNodes(end + 1);
    for (size_t i = 0; i < kdNodes.size(); ++i) {
      kdNodes[i] = new KdNode<K>();
    }
#else
    // Allocate a vector of KdNode instances but only
    // for references that have not been removed.
    tree->kdNodes = new vector<KdNode<K>>(end + 1);
#endif
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    allocateTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

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

    // Fill the permutation vector by calculating the permutation of the indices vector
    // and the the partition coordinate of the tuple at each depth in the tree.
    for (size_t i = 0; i < permutation.size(); ++i) {
      // The last entry of the indices vector contains the partition coordinate.
      indices[numDimensions + 1] = i % numDimensions;
      // Swap the first and second to the last elements of the indices vector.
      swap(indices, 0, numDimensions);
      // Copy the indices vector to one row of the permutation vector.
      permutation[i] = indices;
      // Swap the third and second to the last elements of the indices vector.
      swap(indices, numDimensions - 1, numDimensions);
    }

    // Build the k-d tree with multiple threads if possible.
#ifndef PREALLOCATE
    tree->root = buildKdTree(references, permutation, kdNodes, 0, end, maximumSubmitDepth, 0);
#else
    tree->root = buildKdTree(references, permutation, tree->kdNodes, 0, end, maximumSubmitDepth, 0);
#endif
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    kdTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Verify the k-d tree and report the number of kdNodes. Begin by
    // creating a 1D permutation vector for use by the verifyKdTree function.
    //
    // Because the partition coordinate permutes in the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // (for e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference array
    // using 0 as the leading key of the super key.
    beginTime = steady_clock::now();
    vector<signed_size_t> permutationVerify;
    KdNode<K>::createPermutation(permutationVerify, numDimensions, coordinates.size());
    numberOfNodes = tree->root->verifyKdTree(permutationVerify, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    verifyTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
  
    // Delete the references arrays but not the tuples arrays that they point to.
    beginTime = steady_clock::now();
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // If PREALLOCATE is undefined, clear the kdNodes vector in order to measure its
    // deallocation time. Do not deallocate the KdNode instances to which it points
    // because they will be deallocated recursively by the ~KdNode destructor.
#ifndef PREALLOCATE
    kdNodes.clear();
#endif
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    deallocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Return the pointer to the root of the k-d tree.
    return tree;
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list<KdNode<K>*> that is passed by reference and modified
   * queryLower - the query lower bound vector that is passed by reference and modified
   * queryUpper - the query upper bound vector that is passed by reference and modified
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * size - the number of points in the coordinates vector
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    signed_size_t const size) {

    root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth, size);
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
   * size - the number of points in the coordinates vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    signed_size_t const size,
                    vector<bool> const& enable) {
    
    root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth, size, enable);
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
  void bruteRegion(list<KdNode<K>*>& result,
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
   * size - the number of points in the coordinates vector
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                            signed_size_t const size) {
    
    root->findNearestNeighbors(neighbors, query, numNeighbors, size);
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * size - the number of points in the coordinates vector
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                            signed_size_t const size,
                            vector<bool> const& enable) {
    
    root->findNearestNeighbors(neighbors, query, numNeighbors, size, enable);
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified
   * query - the query vector
   * permutation - vector that specifies permutation of the partition coordinate
   * numNeighbors - the number M of nearest neighbors to attempt to find
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            vector<signed_size_t> const& permutation,
                            signed_size_t const numNeighbors) {
    
    root->findNearestNeighbors(neighbors, query, permutation, numNeighbors);
  }

  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified
   * query - the query vector
   * permutation - vector that specifies permutation of the partition coordinate
   * numNeighbors - the number M of nearest neighbors to attempt to  find
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            vector<signed_size_t> const& permutation,
                            signed_size_t const numNeighbors,
                            vector<bool> const& enable) {

    root->findNearestNeighbors(neighbors, query, permutation, numNeighbors, enable);
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
  void verifyNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighborsFast,
                              forward_list< pair<double, KdNode<K>*> >& neighborsSlow) const {
    
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
   * returns: a sorted forward list of (double, KdNode*) pairs
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  forward_list<pair<double, KdNode<K>*>> sortByDistance(list<KdNode<K>*> const& kdList,
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
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
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
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> >* >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> >* >& rnn,
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
   * numberOfNodes - the number of nodes in the k-d tree
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vectors. Hence, this function is not static.
   */
public:
  void verifyReverseNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                              vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                              signed_size_t const numberOfNodes) {

    root->verifyReverseNeighbors(nn, rnn, numberOfNodes);
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
  void calculateMeanStd(vector< forward_list< pair<double, KdNode<K>*> > >& vec,
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
  size_t nonEmptyLists(vector< forward_list< pair<double, KdNode<K>*> > >& vec) const {

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
  void bruteNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                             vector<K> const& query,
                             signed_size_t const numNeighbors) {
    
    root->bruteNearestNeighbors(neighbors, query, numNeighbors);
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
   * regionList - a list of (double, KdNode*) pairs
   * maximumNumberOfNodesToPrint - the maximum number of KdNodes to print
   * numDimensions - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuples(forward_list<pair<double, KdNode<K>*>> const& regionList,
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
}; // class KdTree

#endif // KD_TREE_KNLOGN_H
