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
 * -D DUAL_THREAD_YUCAO - Execute the dual-threaded version of the k-d tree-building
 *                        algorithm depicted in Figure 2 of Yu Cao et al.
 * 
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * 
 * -D MERGE_CUTOFF=n - A cutoff for switching from 1 to 2 threads to merge reference
 *                     arrays in the KdNode::mergeSort* functions (default 4096)
 * 
 * -D REVERSE_NEAREST_NEIGHBORS - Enable the construction of a reverse nearest neighbors
 *                                list in response to the -r command-line option.
 */

#ifndef KD_TREE_YUCAO_H
#define KD_TREE_YUCAO_H

#include "kdTreeNode.h"

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
   * The buildKdTree function builds a k-d tree via O(n) complexity by recursively subdividing
   * the "final" (aka f) vector via the algorithm described by Yu Cao et al. and depicted in
   * Figure 3 of "A New Method to Construct the KD Tree Based on Presorted Results", Complexity,
   * Volume 2020, Article ID 8883945, https://doi.org/10.1155/2020/8883945
   *
   * Calling parameters:
   *
   * f - a vector<size_t> that represents the "final" vector
   * t - a K** of pointers to each of the (x, y, z, w...) tuples
   * kdNodes - a vector<KdNode<K>>* if PREALLOCATE is defined
   *         - a vector<KdNode<K>*> if PREALLOCATE is undefined
   * size_t init - initial position of the current sub-tree
   * size_t size - size of the current sub-tree
   *
   * returns: a KdNode pointer to the root of the current sub-tree
   */
private:
  static KdNode<K>* buildKdTree(vector<size_t> const& f,
                                K** const t,
#ifndef PREALLOCATE
                                vector<KdNode<K>*> const& kdNodes,
#else
                                vector<KdNode<K>>* const kdNodes,
#endif
                                signed_size_t const init,
                                signed_size_t const size) {

    // Subdivide the current sub-tree according to the number
    // of elements that it contains. The optimum order of the
    // following tests might be determined by counting the
    // number of times that each test is performed for a typical
    // set of tuples.
    KdNode<K>* root = nullptr;
    if (size > 4) {
      // The current sub-tree contains greater than four elements,
      // so its lower and upper sub-trees are created via recursive
      // calls to this buildKdTree function.
      size_t const lower = init;
      size_t const lowerSize = size/2;
      size_t const upperSize = (size - 1)/2;
      size_t const middle = init + lowerSize;
      size_t const upper = 1 + middle;
      root = KdNode<K>::getKdNode(t, kdNodes, f[middle]);
      root->ltChild = buildKdTree(f, t, kdNodes, lower, lowerSize);
      root->gtChild = buildKdTree(f, t, kdNodes, upper, upperSize);
    } else if (size == 4) {
      // The current sub-tree contains four elements, so the element that
      // has the largest key is the upper child of the element that has
      // the next largest key, and the lower child of that element that has
      // the next largest key is a sub-tree that contains the two elements
      // that have the smallest keys.
      size_t const lower = init;
      size_t const lowerSize = size/2;
      size_t const middle = init + lowerSize;
      size_t const upper = 1 + middle;
      root = KdNode<K>::getKdNode(t, kdNodes, f[middle]);
      root->ltChild = buildKdTree(f, t, kdNodes, lower, lowerSize);
      root->gtChild = KdNode<K>::getKdNode(t, kdNodes, f[upper]);
    } else if (size == 3) {
      // The current sub-tree contains three elements, so the middle
      // element is the parent of two 1-element sub-trees.
      size_t const lower = init;
      size_t const middle = init + (size/2);
      size_t const upper = 1 + middle;
      root  = KdNode<K>::getKdNode(t, kdNodes, f[middle]);
      root->ltChild = KdNode<K>::getKdNode(t, kdNodes, f[lower]);
      root->gtChild = KdNode<K>::getKdNode(t, kdNodes, f[upper]);
    } else if (size == 2) {
      // The current sub-tree contains two elements, so the element that
      // has the smaller key is the lower child of the element that
      // has the larger key, and the element that has the larger key
      // has no upper child.
      size_t const lower = init;
      size_t const middle = init + (size/2);
      root = KdNode<K>::getKdNode(t, kdNodes, f[middle]);
      root->ltChild = KdNode<K>::getKdNode(t, kdNodes, f[lower]);
      root->gtChild = nullptr;
    } else if (size == 1) {
      // The current sub-tree comprises only one element,
      // so get a k-d node for it. This particular case
      // might not be executed because the 2, 3, and 4
      // element cases create a 1-element sub-tree. If
      // the 3-element case were removed, this 1-element
      // case would be executed.
      root = KdNode<K>::getKdNode(t, kdNodes, f[init]);
    } else {
      // The current sub-tree contains no elements.
      throw runtime_error("\n\nsub-tree contains no elements\n");
    }

    return root;
  }

  /*
   * The algorithmYuCao1 executes via a single thread the algorithm
   * in Figure 2 of Yu Cao et al.
   * 
   * Calling parameters:
   * 
   * s - a vector<vector<size_t>> of Presorted Results
   * f - a vector<size_t> to return the result of the algorithm
   * n - the number of tuples
   * dim - the number of dimensions (aka k)
   */
private:
  static void algorithmYuCao1(vector<vector<size_t>> const&s,
                              vector<size_t>& f,
                              size_t const n,
                              size_t const dim) {

    // Initialize the bn, ss, and cur vectors as per Figure 1 of Yu Cao et al.
    vector<size_t> bn(n, 0), ss(n, n), cur(n, 0);

#ifdef DEBUG_PRINT
    cout << endl << "bn, ss, and cur vectors initially" << endl;
    for (size_t k = 0; k < n; ++k) {
      cout << bn[k] << "  " << ss[k] << "  " << cur[k] << endl;
    }
    cout << endl;
#endif

    // Execute the algorithm shown in Figure 2 of Yu Cao et al.
    // zeroCount holds the count of zero ss vector elements
    // and signals the end of iterative splitting when it equals
    // the number of tuples.
    for (size_t j = 0, zeroCount = 0; zeroCount < n; ++j) {

      // The reference arrays index advances cyclically with each split.
      size_t d = j % dim;

      // Process each element of the s vector in sorted order.
      for (size_t i = 0; i < n; ++i) {
        size_t tmpI = s[d][i];
        size_t tmpBN = bn[tmpI];
        size_t tmpSize = ss[tmpBN];

        // Don't re-process this element if it is the median of
        // a higher-level sub-tree and hence has already been
        // assigned a size of zero.
        if (tmpSize == 0) {
          continue;
        }

        // Compute the start indices and half-sizes of the lower and
        // upper sub-trees and the median of the current sub-tree,
        // which is in fact the half size of the lower sub-tree.
        // If the current sub-tree contains an odd number of elements,
        // its lower and upper sub-trees contain equal numbers of
        // elements, but if the current sub-tree contains an even
        // number of elements, its upper sub-tree contains one less
        // element than its lower sub-tree.
        size_t lower = tmpBN;
        size_t lowerSize = tmpSize/2;
        size_t upperSize = (tmpSize - 1)/2;
        size_t median = tmpBN + lowerSize;
        size_t upper = 1 + median;

        // Does the element belong to the lower sub-tree?
        if (cur[tmpBN] < lowerSize) {
          bn[tmpI] = lower; // Yes, so assign lower to its beginning.
        }

        // Otherwise, does the element belong to the upper sub-tree?
        else if (cur[tmpBN] > lowerSize) {
          bn[tmpI] = upper; // Yes, so assign upper to its beginning.
        }

        // Otherwise, the element is the median of the current sub-tree,
        // so assign median to its beginning, assign zero to its size, and
        // increment the count of zero ss vector elements.
        else {
          bn[tmpI] = median;
          ss[median] = 0;
          zeroCount++;
        }

        // If processing the current sub-tree has finished, initialize
        // the start indices and sizes of the lower and upper sub-trees
        // unless the size of the current sub-tree is 1 to specify a leaf.
        // No need to initialize cur[tmpBN] because tmpBN equals lower.
        // Presumably, the case of a 2-element sub-tree is accommodated
        // by the padding of the ss and cur vectors by one extra element.
        if (tmpSize > 1 && ++cur[tmpBN] >= tmpSize) {
          ss[lower] = lowerSize;
          ss[upper] = upperSize;
          cur[lower] = cur[upper] = 0;
        }
      }

      // Print the contents of the bn, ss, and cur vectors after each split.
#ifdef DEBUG_PRINT
      cout << endl << "bn, ss, and cur vectors after split " << (j+1) << endl;
      for (size_t k = 0; k < n; ++k) {
        cout << bn[k] << "  " << ss[k] << "  " << cur[k] << endl;
      }
#endif

    }

    // Construct the final (aka f) vector.
    for (size_t i = 0; i < n; ++i) {
      f[bn[i]] = i;
    }

  }

  /*
   * The split function processes each element of the upper half
   * of the s vector in reverse sorted order.
   * 
   * Calling parameters:
   * 
   * s - a vector<vector<size_t> of Presorted Results (aka index vectors)
   * bn - a vector<size_t> of the beginning indices of sub-trees
   * ss - a vector<size_t> of the sizes of sub-trees
   * curFwd - a vector<size_t> of current counts of sub-tree elements
   *          counted in forward manner (low to high address)
   * curRev - a vector<size_t> of current counts of sub-tree elements
   *          counted in reverse manner (high to low address)
   * zeroCount - an atomic<size_t> count of ss zero elements
   * d - the dimension for indexing the Presorted Results (s) vector
   * n - the number of tuples
   */
private:
  static void split(vector<vector<size_t>> const& s,
                    vector<size_t>& bn,
                    vector<size_t>& ss,
                    vector<size_t>& curFwd,
                    vector<size_t>& curRev,
                    atomic<size_t>& zeroCount,
                    size_t const& d,
                    size_t const& n) {

    // Process each element of the s vector in sorted order.
    for (signed_size_t i = n-1; i >= static_cast<signed_size_t>(n/2); --i) {
      size_t tmpI = s[d][i];
      size_t tmpBN = bn[tmpI];
      size_t tmpSize = ss[tmpBN];

      // Don't re-process this element if it is the median of
      // a higher-level sub-tree and hence has already been
      // assigned a size of zero.
      if (tmpSize == 0) {
        continue;
      }

      // Compute the start indices and half-sizes of the lower and
      // upper sub-trees and the median of the current sub-tree,
      // which is in fact the half size of the lower sub-tree.
      // If the current sub-tree contains an odd number of elements,
      // its lower and upper sub-trees contain equal numbers of
      // elements, but if the current sub-tree contains an even
      // number of elements, its upper sub-tree contains one less
      // element than its lower sub-tree.
      size_t lower = tmpBN;
      size_t lowerSize = tmpSize/2;
      size_t upperSize = (tmpSize - 1)/2;
      size_t median = tmpBN + lowerSize;
      size_t upper = 1 + median;

      // Does the element belong to the lower sub-tree?
      if (tmpSize - 1 - curRev[tmpBN] < lowerSize) {
        bn[tmpI] = lower; // Yes, so assign lower to its beginning.
      }

      // Otherwise, does the element belong to the upper sub-tree?
      else if (tmpSize - 1 - curRev[tmpBN] > lowerSize) {
        bn[tmpI] = upper; // Yes, so assign upper to its beginning.
      }

      // Otherwise, the element is the median of the current sub-tree,
      // so assign median to its beginning, assign zero to its size, and
      // increment the count of zero ss vector elements. Apparently, the
      // following fetch_add function is the fastest way to implement
      // zeroCount++ but std::memory_order_relaxed doesn't improve the
      // speed on Intel architectures. See answer #10 at the following URL.
      // https://stackoverflow.com/questions/41206861/how-to-implement-an-atomic-counter
      else {
        bn[tmpI] = median;
        ss[median] = 0;
        zeroCount.fetch_add(1, memory_order_relaxed);
      }

      // If processing the current sub-tree has finished, initialize
      // the start indices and sizes of the lower and upper sub-trees
      // unless the size of the current sub-tree is 1 to specify a leaf.
      // No need to initialize cur[tmpBN] because tmpBN equals lower.
      // Presumably, the case of a 2-element sub-tree is accommodated
      // by the padding of the ss and cur vectors by one extra element.
      if (tmpSize > 1 && ++curRev[tmpBN] + curFwd[tmpBN] >= tmpSize) {
        ss[lower] = lowerSize;
        ss[upper] = upperSize;
        curRev[lower] = curRev[upper] = curFwd[lower] = curFwd[upper] = 0;
      }
    }
  }

  /*
   * The algorithmYuCao2 executes via two threads the algorithm
   * in Figure 2 of Yu Cao et al.
   * 
   * Calling parameters:
   * 
   * s - a vector<vector<size_t>> of Presorted Results
   * f - a vector<size_t> to return the result of the algorithm
   * n - the number of tuples
   * dim - the number of dimensions (aka k)
   */
private:
  static void algorithmYuCao2(vector<vector<size_t>> const&s,
                              vector<size_t>& f,
                              size_t const n,
                              size_t const dim) {

    // Initialize the bn, ss, curFwd, and curRev vectors as per
    // Figure 1 of Yu Cao et al.
    vector<size_t> bn(n, 0), ss(n, n), curFwd(n, 0), curRev(n, 0);

    // Execute the algorithm shown in Figure 2 of Yu Cao et al.
    // zeroCount holds the count of zero ss vector elements
    // and signals the end of iterative splitting when it equals
    // the number of tuples.
    atomic<size_t> zeroCount{0};
    for (size_t j = 0; zeroCount < n; ++j) {

      // The Presorted Results (s) dimension advances cyclically with each split.
      size_t d = j % dim;

      // Process each element of the lower half of the s vector
      // in reverse sorted order via a child thread. This async
      // function call won't compile without a lambda function.
      auto splitFuture = async(launch::async, [&] {
                                                    split(s, bn, ss, curFwd, curRev,
                                                          zeroCount, d, n);
                                                  });

      // Process each element of the lower half of the s vector
      // in forward sorted order via the master thread.
      for (size_t i = 0; i < n/2; ++i) {
        size_t tmpI = s[d][i];
        size_t tmpBN = bn[tmpI];
        size_t tmpSize = ss[tmpBN];

        // Don't re-process this element if it is the median of
        // a higher-level sub-tree and hence has already been
        // assigned a size of zero.
        if (tmpSize == 0) {
          continue;
        }

        // Compute the start indices and half-sizes of the lower and
        // upper sub-trees and the median of the current sub-tree,
        // which is in fact the half size of the lower sub-tree.
        // If the current sub-tree contains an odd number of elements,
        // its lower and upper sub-trees contain equal numbers of
        // elements, but if the current sub-tree contains an even
        // number of elements, its upper sub-tree contains one less
        // element than its lower sub-tree.
        size_t lower = tmpBN;
        size_t lowerSize = tmpSize/2;
        size_t upperSize = (tmpSize - 1)/2;
        size_t median = tmpBN + lowerSize;
        size_t upper = 1 + median;

        // Does the element belong to the lower sub-tree?
        if (curFwd[tmpBN] < lowerSize) {
          bn[tmpI] = lower; // Yes, so assign lower to its beginning.
        }

        // Otherwise, does the element belong to the upper sub-tree?
        else if (curFwd[tmpBN] > lowerSize) {
          bn[tmpI] = upper; // Yes, so assign upper to its beginning.
        }

        // Otherwise, the element is the median of the current sub-tree,
        // so assign median to its beginning, assign zero to its size, and
        // increment the count of zero ss vector elements. Apparently, the
        // following fetch_add function is the fastest way to implement
        // zeroCount++ but std::memory_order_relaxed doesn't improve the
        // speed on Intel architectures. See answer #10 at the following URL.
        // https://stackoverflow.com/questions/41206861/how-to-implement-an-atomic-counter
        else {
          bn[tmpI] = median;
          ss[median] = 0;
          zeroCount.fetch_add(1, memory_order_relaxed);
        }

        // If processing the current sub-tree has finished, initialize
        // the start indices and sizes of the lower and upper sub-trees
        // unless the size of the current sub-tree is 1 to specify a leaf.
        // No need to initialize cur[tmpBN] because tmpBN equals lower.
        // Presumably, the case of a 2-element sub-tree is accommodated
        // by the padding of the ss and cur vectors by one extra element.
        if (tmpSize > 1 && ++curFwd[tmpBN] + curRev[tmpBN] >= tmpSize) {
          ss[lower] = lowerSize;
          ss[upper] = upperSize;
          curFwd[lower] = curFwd[upper] = curRev[lower] = curRev[upper] = 0;
        }
      }

      // Wait for the child thread to finish execution.
      try {
        splitFuture.get();
      }
      catch (exception const& e) {
        throw runtime_error("\n\ncaught exception for split future in algrithmYuCao2\n");
      }

      // Print the contents of the bn, ss, curFWd, and curRev vectors after each split.
#ifdef DEBUG_PRINT
      cout << endl << "bn, ss, curFwd, and curRev vectors after split " << (j+1) << endl;
      for (size_t k = 0; k < n; ++k) {
        cout << bn[k] << "  " << ss[k] << "  " << curFwd[k] << "  " << curRev[k] << endl;
      }
#endif

    }

    // Construct the final (aka f) vector.
    for (size_t i = 0; i < n; ++i) {
      f[bn[i]] = i;
    }
  }

  /*
   * The createKdTree function presorts the references arrays via O(kn log n) complexity
   * and then conceptually splits those reference arrays via O(n log n) complexity to create
   * the "final" (aka f) vector via the algorithm described by Yu Cao et al. and depicted in
   * Figure 2 of "A New Method to Construct the KD Tree Based on Presorted Results", Complexity,
   * Volume 2020, Article ID 8883945, https://doi.org/10.1155/2020/8883945
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

    // Allocate the reference arrays including one additional array.
    auto beginTime = steady_clock::now();
    size_t n = coordinates.size();
    size_t numDimensions = coordinates[0].size();
    K*** references = new K**[numDimensions + 1];
    for (size_t i = 0; i < numDimensions + 1; ++i) {
      references[i] = new K*[n];
    }

#ifndef PREALLOCATE
    // Allocate tuple arrays for the first reference array.
    // Each tuple array will be deallocated by either the
    // KdNode::removeDuplicates function or the ~KdNode destructor.
    // Copy (x, y, z, w...) coordinates to each tuple array.
    // Create an additional element in each tuple array to
    // store the index of the tuple as required by Yu Cao's
    // algorithm. This kludge will work for an integer tuple
    // array but perhaps not for a floating-point array.
    for (size_t i = 0; i < n; ++i) {
      references[0][i] = new K[numDimensions+1];
      for (size_t j = 0; j < numDimensions; ++j) {
        references[0][i][j] = coordinates[i][j];
      }
      references[0][i][numDimensions] = i;
    }

#ifdef DEBUG_PRINT
    cout << "references[0] array" << endl;
    for (size_t i = 0; i < n; ++i) {
      cout << "(" << references[0][i][0] << "," << references[0][i][1] << ") " << references[0][i][2] << endl;
    }
    cout << endl;
#endif

#else
    // Allocate a tuple vector that comprises all tuple arrays and
    // assign individual tuple arrays to the first reference array.
    // This vector will be deallocated by the ~KdTree destructor.
    // Copy (x, y, z, w...) coordinates to each tuple array.
    // Create an additional element in each tuple array to
    // store the index of the tuple as required by Yu Cao's
    // algorithm. This kludge will work for an integer tuple
    // array but perhaps not for a floating-point array.
    tree->tuples = new vector<K>((numDimensions+1) * n);
    for (size_t i = 0; i < n; ++i) {
      references[0][i] = &(*(tree->tuples)).data()[(numDimensions+1) * i];
      for (size_t j = 0; j < numDimensions; ++j) {
        references[0][i][j] = coordinates[i][j];
      }
      references[0][i][numDimensions] = i;
    }

#ifdef DEBUG_PRINT
    cout << "references[0] array" << endl;
    for (size_t i = 0; i < n; ++i) {
     cout << "(" << references[0][i][0] << "," << references[0][i][1] << ") " << references[0][i][2] << endl;
    }
    cout << endl;
#endif

#endif
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    allocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Sort the first reference array using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
    // Also, only the first reference array is populated with tuple arrays.
    beginTime = steady_clock::now();
    MergeSort<K>::mergeSortReferenceAscending(references[0], references[numDimensions], 0,
                                              n-1, 0, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    sortTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

#ifdef DEBUG_PRINT
    cout << "references[0] array" << endl;
    for (size_t i = 0; i < n; ++i) {
     cout << "(" << references[0][i][0] << "," << references[0][i][1] << ") " << references[0][i][2] << endl;
    }
    cout << endl;

    cout << "coordinates vector" << endl;
    for (size_t i = 0; i < n; ++i) {
      cout << "(" << coordinates[i][0] << "," << coordinates[i][1] << ")" << endl;
    }
    cout << endl;
#endif

    // Remove references to duplicate coordinates via one pass through the first reference array
    // and adjust the number of coordinates (n) to represent the remaining references. 
    beginTime = steady_clock::now();
    signed_size_t const end = KdNode<K>::removeDuplicates(references[0], 0, numDimensions, n);
    n = end + 1;
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    removeTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Verify that no tuple's index lies beyond the range [0, end];
    for (size_t i = 0; i < n; ++i) {
      if (references[0][i][numDimensions] > end) {
        ostringstream buffer;
        buffer << "\n\ntuple " << i << " index = " << references[0][i][numDimensions]
               << " > " << end << endl << endl;
        throw runtime_error(buffer.str());
      }
    }

#ifdef DEBUG_PRINT
    cout << "end = " << end << endl;
    cout << "references[0] array" << endl;
    for (size_t i = 0; i < n; ++i) {
     cout << "(" << references[0][i][0] << "," << references[0][i][1] << ") " << references[0][i][2] << endl;
    }
    cout << endl;
#endif

    // Copy the tuple array pointers from first reference array to all but
    // the last reference array. Sort these referece arrays using multiple
    // threads, but only to the extent of those tuples that were not "removed",
    // i.e., skipped by the removeDuplicates function above.
    beginTime = steady_clock::now();
    for (size_t i = 1; i < numDimensions; ++i) {
      for (size_t j = 0; j < n; ++j) {
        references[i][j] = references[0][j];
      }

#ifdef DEBUG_PRINT
      cout << "references[" << i << "] array before merge sort" << endl;
      for (size_t k = 0; k < n; ++k) {
      cout << "(" << references[i][k][0] << "," << references[i][k][1] << ") " << references[i][k][2] << endl;
      }
      cout << endl;
#endif

      MergeSort<K>::mergeSortReferenceAscending(references[i], references[numDimensions], 0,
                                                n-1, i, numDimensions, maximumSubmitDepth, 0);
    }
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    sortTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

#ifdef DEBUG_PRINT
    cout << "references[1] array after merge sort" << endl;
    for (size_t i = 0; i < n; ++i) {
     cout << "(" << references[1][i][0] << "," << references[1][i][1] << ") " << references[1][i][2] << endl;
    }
    cout << endl;
#endif

    beginTime = steady_clock::now();
#ifndef PREALLOCATE
    // Allocate the kdNodes vector of pointers to KdNode instances
    // but only for references that were not removed by removeDuplicates.
    // Allocate a KdNode instance for each entry in the kdNodes vector.
    vector<KdNode<K>*> kdNodes(n);
    for (size_t i = 0; i < kdNodes.size(); ++i) {
      kdNodes[i] = new KdNode<K>();
    }
#else
    // Allocate a vector of KdNode instances for references that
    // were not removed by the removeDuplicates function,
    tree->kdNodes = new vector<KdNode<K>>(n);
#endif
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    allocateTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Start the timer to time building the k-d tree.
    beginTime = steady_clock::now();

    // Create Presorted Results (aka s) vectors that contain the
    // indices of tuples, as required by Yu Cao's algorithm.
    vector<vector<size_t>> s(numDimensions, vector<size_t>(n));
    for (size_t i = 0; i < numDimensions; ++i) {
      for (size_t j = 0; j < n; ++j) {
        s[i][j] = references[i][j][numDimensions];
      }
    }

#ifdef DEBUG_PRINT
    cout << "s vector " << coordinates.size() << endl;
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < numDimensions; ++j) {
        cout << s[j][i] << "  ";
      }
      cout << endl;
    }
    cout << endl;
#endif

    // Call either the single-threaded or dual-threaded version of the
    // algorithm depicted in Figure 2 of Yu Cao et al., depending on
    // whether 1 or 2 threads are available to execute the algorithm.
    // Create the f vector to return the result of either version.
    vector<size_t> f(n);
#ifdef DUAL_THREAD_YUCAO
    if (maximumSubmitDepth < 0) {
      algorithmYuCao1(s, f, n, numDimensions);
    } else {
      algorithmYuCao2(s, f, n, numDimensions);
    }
#else
    algorithmYuCao1(s, f, n, numDimensions);
#endif

#ifdef DEBUG_PRINT
    cout << endl << "final vector" << endl;
    for (size_t i = 0; i < n; ++i) {
      cout << f[i] << endl;
    }
    cout << endl;
#endif

    // Restore the original, unsorted order of the first reference
    // array into the last reference array, because Yu Cao's algorithm
    // expects unsorted tuples.
    for (size_t i = 0; i < n; ++i) {
      references[numDimensions][references[0][i][numDimensions]] = references[0][i];
    }

#ifdef DEBUG_PRINT
    cout << "references[2] array" << endl;
    for (size_t i = 0; i < n; ++i) {
     cout << "(" << references[2][i][0] << "," << references[2][i][1] << ") " << references[2][i][2] << endl;
    }
    cout << endl;
#endif

    // Build the k-d tree from the "final" (aka f) vector.
#ifndef PREALLOCATE
    tree->root = buildKdTree(f, references[numDimensions], kdNodes, 0, n);
#else
    tree->root = buildKdTree(f, references[numDimensions], tree->kdNodes, 0, n);
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
    for (size_t i = 0; i < numDimensions+1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // If PREALLOCATE is undefined, clear the kdNodes vector in order to measure its
    // deallocation time. Do not deallocate the KdNode instances that it points to
    // because they will be deallocated recursively by the ~KdNode destructor.
#ifndef PREALLOCATE
    kdNodes.clear();
#endif

    // Clear the s and f vectors in order to measure their deallocation times.
    s.clear();
    f.clear();

    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    deallocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

#ifdef DEBUG_PRINT
    tree->printKdTree(numDimensions, 0);
    cout << endl;
#endif

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

#endif // KD_TREE_YUCAO_H
