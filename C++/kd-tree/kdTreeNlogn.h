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
 * -D MEDIAN_OF_MEDIANS_CUTOFF=n - A cutoff for switching from median of medians to
 *                                 insertion sort in KdTree::partition (default 15)
 * 
 * -D MEDIAN_CUTOFF=n - A cutoff for switching from 1 to 2 threads to calculate the median
 *                      in the KdTree::partition function (default 16384)
 * 
 * -D INDEX_CUTOFF=n - A cutoff for switching from 1 to 2 threads to find the index of
 *                     the calculated median in KdTree::partition (default 16384)
 * 
 * -D BIDIRECTIONAL_PARTITION - Partition an array about the median of medians proceeding
 *                              from both ends of the array instead of only the beginning.
 * 
 * -D NLOGN_CUTOFF=n - A cutoff for using multiple threads in buildKdTree (default 4096)
 */

#ifndef KD_TREE_NLOGN_H
#define KD_TREE_NLOGN_H

/* A cutoff for using multiple threads in buildKdTree */
#ifndef NLOGN_CUTOFF
#define NLOGN_CUTOFF (4096)
#endif

/* A cutoff for switching from median of medians to insertion sort in KdNode::partition */
#ifndef MEDIAN_OF_MEDIANS_CUTOFF
#define MEDIAN_OF_MEDIANS_CUTOFF (15)
#endif

/* A cutoff for switching from 1 to 2 threads to calculate the median in KdNode::partition */
#ifndef MEDIAN_CUTOFF
#define MEDIAN_CUTOFF (16384)
#endif

/* A cutoff for switching from 1 to 2 threads to find the index of the median in KdNode::partition */
#ifndef INDEX_CUTOFF
#define INDEX_CUTOFF (1073741824) // =2^30 to disable switching to 2 threads
#endif

#include "kdTree.h"

/* Forward references to all classes to support any order of compilation */
template <typename>
class KdTreeDynamic;

template <typename>
class KdTree;

template <typename>
class KdNode;

template <typename>
class MergeSort;

template <typename>
class NearestNeighborHeap;

template <typename>
class KdTreeNlogn;

template <typename>
class KdTreeKnlogn;

template <typename>
class KdTreeYuCao;

/* The KdTreeNlogn class provides static factory functions. */
template <typename K>
class KdTreeNlogn
{
  /*
   * The swap function swaps two array elements.
   *
   * Calling parameters:
   *
   * a - K** array of references to the (x,y,z,w...) coordinates
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  inline
  static void swap(K** const a,
                   signed_size_t const i,
                   signed_size_t const j) {
    
    auto const t = a[i];
    a[i] = a[j];
    a[j] = t;
  }

  /*
   * The following select_j_k functions select the jth of k items.  Adapted
   * from Chapter 4, "Linear Orderings", of Alexander Stepanov's and
   * Paul McJones' "Elements of Programming", Addison-Wesley, New York, 2009.
   *
   * Calling parameters:
   *
   * a through e - arrays to include in the selection
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   *
   * returns - a K* that represents the selected array
   */
private:
  inline
  static K* select_0_2(K* const a,
                       K* const b,
                       signed_size_t const p,
                       signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(a, b, p, dim) < 0) {
      // a < b
      return a;
    }
    else {
      // b < a
      return b;
    }
  }

private:
  inline
  static K* select_1_2(K* const a,
                       K* const b,
                       signed_size_t const p,
                       signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(a, b, p, dim) < 0) {
      // a < b
      return b;
    }
    else {
      // b < a
      return a;
    }
  }

private:
  inline
  static K* select_1_3_ab(K* const a,
                          K* const b,
                          K* const c,
                          signed_size_t const p,
                          signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(b, c, p, dim) < 0) {
      // a < b < c
      return b;
    }
    else {
      // a ? c < b
      return select_1_2(a, c, p, dim);
    }
  }

private:
  inline
  static K* select_1_3(K* const a,
                       K* const b,
                       K* const c,
                       signed_size_t const p,
                       signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(a, b, p, dim) < 0) {
      // a < b
      return select_1_3_ab(a, b, c, p, dim);
    }
    else {
      // b < a
      return select_1_3_ab(b, a, c, p, dim);
    }
  }

private:
  inline
  static K* select_1_4_ab_cd(K* const a,
                             K* const b,
                             K* const c,
                             K* const d,
                             signed_size_t const p,
                             signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(c, a, p, dim) < 0) {
      // c < a < b && a ? d so c is eliminated and a ? d
      return select_0_2(a, d, p, dim);
    }
    else {
      // a < b ? c < d so a is eliminated and b ? c
      return select_0_2(b, c, p, dim);
    }
  }

private:
  inline
  static K* select_1_4_ab(K* const a,
                          K* const b,
                          K* const c,
                          K* const d,
                          signed_size_t const p,
                          signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(c, d, p, dim) < 0) {
      // a < b && c < d
      return select_1_4_ab_cd(a, b, c, d, p, dim);
    }
    else {
      // a < b && d < c
      return select_1_4_ab_cd(a, b, d, c, p, dim);
    }
  }

private:
  inline
  static K* select_1_4(K* const a,
                       K* const b,
                       K* const c,
                       K* const d,
                       signed_size_t const p,
                       signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(a, b, p, dim) < 0) {
      // a < b
      return select_1_4_ab(a, b, c, d, p, dim);
    }
    else {
      // b < a
      return select_1_4_ab(b, a, c, d, p, dim);
    }
  }

private:
  inline
  static K* select_2_5_ab_cd(K* const a,
                             K* const b,
                             K* const c,
                             K* const d,
                             K* const e,
                             signed_size_t const p,
                             signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(c, a, p, dim) < 0) {
      // c < a < b && c < d ? e so c is eliminated and a < b && d ? e
      return select_1_4_ab(a, b, d, e, p, dim);
    }
    else {
      // a < b ? c && c < d ? e && b ? e so a is eliminated and c < d && b ? e
      return select_1_4_ab(c, d, b, e, p, dim);
    }
  }

private:
  inline
  static K* select_2_5_ab(K* const a,
                          K* const b,
                          K* const c,
                          K* const d,
                          K* const e,
                          signed_size_t const p,
                          signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(c, d, p, dim) < 0) {
      // a < b && c < d
      return select_2_5_ab_cd(a, b, c, d, e, p, dim);
    }
    else {
      // a < b && d < c
      return select_2_5_ab_cd(a, b, d, c, e, p, dim);
    }
  }

private:
  inline
  static K* select_2_5(K* const a,
                       K* const b,
                       K* const c,
                       K* const d,
                       K* const e,
                       signed_size_t const p,
                       signed_size_t const dim) {
    
    if (MergeSort<K>::superKeyCompare(a, b, p, dim) < 0) {
      // a < b
      return select_2_5_ab(a, b, c, d, e, p, dim);
    }
    else {
      // b < a
      return select_2_5_ab(b, a, c, d, e, p, dim);
    }
  }

  /*
   * The partition function partitions an array of references to (x,y,z,w...)
   * tuples about its kth element and returns the array index of the kth element.
   * It implements the algorithm described by Manuel Blum, et al. in "Time Bounds
   * for Selection" in the Journal of Computer and System Sciences, 7:448-461, 1973.
   *
   * See also https://yiqi2.wordpress.com/2013/07/03/median-of-medians-selection-algorithm/
   * that contains a bug that causes a java.lang.ArrayIndexOutOfBoundsException.
   *
   * a - array of references to the (x,y,z,w...) coordinates
   * start - the start index for the elements to be considered
   * n - the number of elements to consider
   * size - the size of the array of references
   * k - the element to find
   * medians - a scratch array for the medians
   * first - the first index for the scratch array
   * p - the most significant dimension or the partition coordinate
   * dim - the number of dimensions
   * twoThreads - use two threads for median calculation
   *
   * returns - the index of the kth element in the array about which the array has been partitioned
   */
private:
  static signed_size_t partition(K** const a,
                                 signed_size_t const start,
                                 signed_size_t const n,
                                 signed_size_t const size,
                                 signed_size_t const k,
                                 K** const medians,
                                 signed_size_t const first,
                                 signed_size_t const p,
                                 signed_size_t const dim,
                                 bool const twoThreads) {

    if (n <= 0 || n > size) {
      ostringstream buffer;
      buffer << "\n\nn = " << n << "  size = " << size << " in partition\n";
      throw runtime_error(buffer.str());
    }
    if (k <= 0 || k > n) {
      ostringstream buffer;
      buffer << "\n\nk = " << k << " in partition\n";
      throw runtime_error(buffer.str());
    }
    if (start + n > size) {
      ostringstream buffer;
      buffer << "\n\nstart = " << start << "  n = " << n << "  size = " << size << " in partition\n";
    }

    // This trivial case terminates recursion.
    if (n == 1 && k == 1) {
      return start;
    }

    // Use insertion sort instead of the median of medians algorithm for a small number of elements,
    // via Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
    // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the array a.
    if (n <= MEDIAN_OF_MEDIANS_CUTOFF) {
      for (signed_size_t i = start + 1; i <= start + n - 1; ++i) {
        auto const tmp = a[i];
        signed_size_t j;
        for (j = i; j > start && MergeSort<K>::superKeyCompare(a[j - 1], tmp, p, dim) > 0; --j) {
          a[j] = a[j - 1];
        }
        a[j] = tmp;
      }
      return start + k - 1;
    }

    // Otherwise, determine how many medians to find.  Round down to count
    // only groups that comprise fully GROUP_SIZE elements.  Any remaining
    // group of elements that doesn't comprise GROUP_SIZE elements will
    // be processed after the following 'for' loop.
    signed_size_t const GROUP_SIZE = 5; // Must be 5 due to select_2_5 function below.
    signed_size_t m = n / GROUP_SIZE;
    signed_size_t startOfGroup;

    // Is more than one thread available to calculate the medians, and are
    // there sufficient medians to justify multi-threaded processing?
    if (twoThreads && m > MEDIAN_CUTOFF) {

      // Yes, calculate the relative index of the middle median.
      signed_size_t const mid = (m + 1) >> 1;
      startOfGroup = mid * GROUP_SIZE;

      // Calculate the lower set of medians with a child thread.
      auto medianFuture =
        async(launch::async, [&] {
                               for (signed_size_t firstOfGroup = 0, i = 0; i < mid; ++i) {

                                 // Find the median of the group of GROUP_SIZE elements via select_2_5.
                                 medians[first + i] = select_2_5(a[start + firstOfGroup],
                                                                 a[start + firstOfGroup + 1],
                                                                 a[start + firstOfGroup + 2],
                                                                 a[start + firstOfGroup + 3],
                                                                 a[start + firstOfGroup + 4],
                                                                 p,
                                                                 dim);

                                 // Update the index to the next group of GROUP_SIZE elements.
                                 firstOfGroup += GROUP_SIZE;
                               }
                             });

      // Calculate the upper set of medians with the current thread.
      for (signed_size_t i = mid; i < m; ++i) {

        // Find the median of the group of GROUP_SIZE elements via select_2_5.
        medians[first + i] = select_2_5(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        a[start + startOfGroup + 4],
                                        p,
                                        dim);

        // Update the index to the next group of GROUP_SIZE elements.
        startOfGroup += GROUP_SIZE;
      }

      // Wait for the child thread to finish execution.
      try {
        medianFuture.get();
      }
      catch (exception const& e) {
        throw runtime_error("\n\ncaught exception for median future in partition\n");
      }
    }
    else
    {
      // No, only one thread is available, so calculate all medians with the current thread.
      startOfGroup = 0;
      for (signed_size_t i = 0; i < m; ++i) {

        // Find the median of the group of GROUP_SIZE elements via select_2_5.
        medians[first + i] = select_2_5(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        a[start + startOfGroup + 4],
                                        p,
                                        dim);

        // Update the index to the next group of GROUP_SIZE elements.
        startOfGroup += GROUP_SIZE;
      }
    }

    // Calculate and check the number of remaining elements.
    signed_size_t const remainingElements = n - startOfGroup;
    if (remainingElements < 0 || remainingElements >= GROUP_SIZE) {
      throw runtime_error("\n\nincorrect group calculation in partition\n");
    }

    // Find the median of any remaining elements via select_j_k.
    switch (remainingElements) {
      case 0:
        break;
      case 1:
        medians[first + m] = a[start + startOfGroup];
        ++m;
        break;
      case 2:
        medians[first + m] = select_0_2(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        p,
                                        dim);
        ++m;
        break;
      case 3:
        medians[first + m] = select_1_3(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        p,
                                        dim);
        ++m;
        break;
      case 4:
        medians[first + m] = select_1_4(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        p,
                                        dim);
        ++m;
        break;
      default:
      {
        ostringstream buffer;
        buffer << "\n\nunhandled case in switch: remainingElements = " << remainingElements << " in partition\n";
        throw runtime_error(buffer.str());
      }
    }

    // Select the median of medians for partitioning the elements.  Note that (m + 1) >> 1
    // correctly designates the median element as the "kth" element instead of the address
    // of the median element in the medians array.  The medians array must start with element
    // first + m for the next level of recursion to avoid overwriting the median elements.
    // The medians array will have adequate capacity for all of the medians for all of the
    // recursive calls because the initial call to this partition function from the buildKdTree
    // function provides a medians array that is the same size as the reference array.  Each
    // recursive call creates the (1 / GROUP_SIZE) fraction of the medians as the call at the
    // prior level of recursion, so the total requirement for storage of medians is the following
    // fraction of the temporary array for the following values of GROUP_SIZE:
    //
    // for GROUP_SIZE = 3, the fraction is 1/3 + 1/9 + 1/27 + 1/81 + ... < 1/2
    // for GROUP_SIZE = 5, the fraction is 1/5 + 1/25 + 1/125 + 1/625 + ... < 1/4
    // for GROUP_SIZE = 7, the fraction is 1/7 + 1/49 + 1/343 + 1/2,401 + ... < 1/6
    // for GROUP_SIZE = 9, the fraction is 1/9 + 1/81 + 1/729 + 1/6,561 + ... < 1/8
    //
    // Note: it is possible to allocate the medians array locally to this partition method
    // instead of providing it via a calling parameter to this method; however, because the
    // mergeSort method requires a temporary array, that array is re-used as the medians array.
    auto const* const medianOfMedians =
      medians[partition(medians, first, m, first + m, (m + 1) >> 1, medians, first + m, p, dim, twoThreads)];

    // Find the index of the median of medians and swap it into a[start + n - 1]
    // so that it is not examined during partitioning of the array a.
    //
    // Is more than one thread available to find the index of the median of medians
    // and are there sufficient array elements to justify dual-threaded processing?
    if (twoThreads && n > INDEX_CUTOFF)
    {
      // Yes, more than one thread is available, so calculate the relative index of the middle element.
      signed_size_t const middle = (n + 1) >> 1;

      // Search for the index in the lower half of the array a with a child thread.
      auto indexFuture =
        async(launch::async, [&] {
                               for (signed_size_t i = 0; i < middle; ++i) {
                                 if (a[start + i] == medianOfMedians) {
                                   swap(a, start + i, start + n - 1);
                                   break;
                                 }
                               }
                             });

      // Search for the index in the upper half of the array a with the current thread.
      for (signed_size_t i = middle; i < n - 1; ++i) {
        if (a[start + i] == medianOfMedians) {
          swap(a, start + i, start + n - 1);
          break;
        }
      }

      // Wait for the child thread to finish execution.
      try {
        indexFuture.get();
      }
      catch (exception const& e) {
        throw runtime_error("\n\ncaught exception for index future in partition\n");
      }
    }
    else
    {
      // No, there are insufficient array elements to justify dual-threadedd processing,
      // so use only one thread to find the index of the median of medians.
      for (signed_size_t i = 0; i < n - 1; ++i) {
        if (a[start + i] == medianOfMedians) {
          swap(a, start + i, start + n - 1);
          break;
        }
      }
    }

    // Partition the array a relative to the median of medians into < and > subsets.
    signed_size_t i = 0;

#ifdef BIDIRECTIONAL_PARTITION
    // Search from both ends of the array a in order to minimize the use of the swap
    // function at the expense of greater use of the MergeSort::superKeyCompare function.
    //
    // NOTE, however, that searching from both ends appears to degrade performance.
    signed_size_t j = n - 2;
    while (j > i) {
      if (MergeSort<K>::superKeyCompare(a[start + i], medianOfMedians, p, dim) < 0) {
        ++i;
      }
      else if (MergeSort<K>::superKeyCompare(a[start + j], medianOfMedians, p, dim) > 0) {
        --j;
      }
      else {
        swap(a, start + i, start + j);
        ++i;
        --j;
      }
    }

    // Ensure that all elements of the < subset are located below a[start + i].
    for (; i < n - 1; ++i) {
      if (MergeSort<K>::superKeyCompare(a[start + i], medianOfMedians, p, dim) > 0) {
        break;
      }
    }
#else
    // Search upward from the beginning of the array a in order to minimize the use of
    // the MergeSort<K>::superKeyCompare function at the expense of greater use of the swap function.
    for (signed_size_t j = 0; j < n - 1; ++j) {
      if (MergeSort<K>::superKeyCompare(a[start + j], medianOfMedians, p, dim) < 0) {
        if (j != i) {
          swap(a, start + j, start + i);
        }
        ++i;
      }
    }
#endif

    // Swap the median of medians into a[start + i] between the < and > subsets.
    swap(a, start + i, start + n - 1);

    // k is 1-based but i is 0-based, so compare k to i + 1 and
    // determine which subset (if any) must be partitioned recursively.
    if (k < i + 1) {

      // The median of medians occupies a position below i, so partition
      // the array elements of the < subset; for this subset, the
      // original kth element is still the kth element of this subset.
      return partition(a, start, i, size, k, medians, first, p, dim, twoThreads);

    }
    else if (k > i + 1) {

      // The median of medians occupies a position above i, so partition
      // the array elements of the > subset; for this subset, the
      // original kth element is not the kth element of this subset
      // because i + 1 elements are in the < subset.
      return partition(a, start + i + 1, n - i - 1, size, k - i - 1,
                      medians, first, p, dim, twoThreads);

    }
    else {

      // The median of medians occupies a[start + i] because k == i + 1, so no
      // further partitioning is necessary.  Return start + i as the index of
      // the kth element under the definition that start is the zeroth index.
      return start + i;
    }
  }

  /*
   * The  buildKdTree function builds a k-d tree by recursively partitioning the reference
   * array and adding KdNodes to the tree.  The super key to be used for partitioning is permuted
   * cyclically for successive levels of the tree in order that sorting use x, y, z, etc. as the
   * most significant portion of the super key.
   *
   * Calling parameters:
   *
   * reference - a K** of references to recursively partition via its K* (x, y, z, w...) tuples
   * temporary - a K** scratch array into which to copy references
   * permutation - a vector<signed_size_t> that indications permutation of the partition coordinate
   * kdNodes - a vector<KdNode<K>>* if PREALLOCATE is defined and KD_TREE_DYNAMIC_h is undefined
   *         - a vector<KdNode<K>*> if PREALLOCATE is undefined or KD_TREE_DYNAMIC_H is defined
   * start - start element of the reference array
   * end - end element of the reference array
   * size - the size of the reference array
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K>* buildKdTree(K** const reference,
                                K** const temporary,
                                vector<signed_size_t> const& permutation,
#if !defined(PREALLOCATE) || defined(KD_TREE_DYNAMIC_H)
                                vector<KdNode<K>*> const& kdNodes,
#else
                                vector<KdNode<K>>* const kdNodes,
#endif
                                signed_size_t const start,
                                signed_size_t const end,
                                signed_size_t const size,
                                signed_size_t const dim,
                                signed_size_t const maximumSubmitDepth,
                                signed_size_t const depth) {

    KdNode<K>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies leading dimension of the key.
    signed_size_t const p = permutation[depth];

    if (end == start) {

      // Only one reference was passed to this method, so store it at this level of the tree.
      node = KdNode<K>::getKdNode(reference, kdNodes, start);
#ifdef KD_TREE_DYNAMIC_H
      node->height = 1;
#endif

    }
    else if (end == start + 1) {

      // Two references were passed to this method in unsorted order, so store the
      // start reference at this level of the tree and determine whether to store the
      // end reference as the < child or the > child.
      node = KdNode<K>::getKdNode(reference, kdNodes, start);
      if (MergeSort<K>::superKeyCompare(reference[start], reference[end], p, dim) > 0) {
        node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, end);
#ifdef KD_TREE_DYNAMIC_H
        node->ltChild->height = 1;
        node->height = 2;
#endif
      }
      else {
        node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);
#ifdef KD_TREE_DYNAMIC_H
        node->gtChild->height = 1;
        node->height = 2;
#endif
      }

    }
    else if (end == start + 2) {

      // Three references were passed to this method in unsorted order, so compare
      // the three references to determine which reference is the median reference.
      // Store the median reference at this level of the tree, store the smallest
      // reference as the < child and store the largest reference as the > child.
      signed_size_t mid = start + 1;
      if (MergeSort<K>::superKeyCompare(reference[start], reference[mid], p, dim) < 0) {
        // reference[start] < reference[mid]
        if (MergeSort<K>::superKeyCompare(reference[mid], reference[end], p, dim) < 0) {
          // reference[start] < reference[mid] < reference[end]
          node = KdNode<K>::getKdNode(reference, kdNodes, mid);
          node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, start);
          node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);
        }
        else {
          // reference[start] < reference[mid]; reference[end] < reference[mid]
          if (MergeSort<K>::superKeyCompare(reference[start], reference[end], p, dim) < 0) {
            // reference[start] < reference[end] < reference[mid]
            node = KdNode<K>::getKdNode(reference, kdNodes, end);
            node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, start);
            node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, mid);
          }
          else {
            // reference[end] < reference[start] < reference[mid]
            node = KdNode<K>::getKdNode(reference, kdNodes, start);
            node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, end);
            node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, mid);
          }
        }
      }
      else {
        // reference[mid] < reference[start]
        if (MergeSort<K>::superKeyCompare(reference[start], reference[end], p, dim) < 0) {
          // reference[mid] < reference[start] < reference[end]
          node = KdNode<K>::getKdNode(reference, kdNodes, start);
          node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, mid);
          node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);
        }
        else {
          // reference[mid] < reference[start]; reference[end] < reference[start]
          if (MergeSort<K>::superKeyCompare(reference[mid], reference[end], p, dim) < 0) {
            // reference[mid] < reference[end] < reference[start]
            node = KdNode<K>::getKdNode(reference, kdNodes, end);
            node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, mid);
            node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, start);
          }
          else {
            // reference[end] < reference[mid] < reference[start]
            node = KdNode<K>::getKdNode(reference, kdNodes, mid);
            node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, end);
            node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, start);
          }
        }
      }
#ifdef KD_TREE_DYNAMIC_H
      node->ltChild->height = node->gtChild->height = 1;
      node->height = 2;
#endif

    }
    else if (end > start + 2) {

      // Four or more references were passed to this method, so calculate the offset
      // of the median element. partition the reference array about its median element,
      // which is the kth element as calculated below.  Note that a < child and
      // a > child will always exist for this case.
      signed_size_t const n = end - start + 1;
      signed_size_t const k = (n + 1) >> 1;

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree, and are
      // there sufficient array elements to justify spawning a child thread?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth
          || end - start < NLOGN_CUTOFF) {

        // No, child threads are not available, so find the median element then
        // partition the reference array about it.  Store the median element
        // from the reference array in a new KdNode.
        signed_size_t const median = partition(reference, start, n, size, k, temporary, start, p, dim, false);
        node = KdNode<K>::getKdNode(reference, kdNodes, median);

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, kdNodes, start,
                                    median - 1, size, dim, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so find the median element then partition
        // the reference array about it.  Store the median element from the reference
        // array in a new KdNode.
        signed_size_t const median = partition(reference, start, n, size, k, temporary, start, p, dim, true);
        node = KdNode<K>::getKdNode(reference, kdNodes, median);

        // Recursively build the < branch of the tree with a child thread.
        auto buildFuture = async(launch::async,
                                buildKdTree,
                                reference,
                                temporary,
                                ref(permutation),
                                ref(kdNodes),
                                start,
                                median - 1,
                                size,
                                dim,
                                maximumSubmitDepth,
                                depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for build future in buildKdTree\n");
        }
      }
#ifdef KD_TREE_DYNAMIC_H
      // Compute the height at this node as the recursion unwinds.
      node->height = KdTreeDynamic<K>::computeHeight(node);
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
   * The buildKdTreePresorted function method builds a k-d tree by using
   * the median of the pre-sorted reference array to partition that array,
   * then calls the buildKdTree function to recursively partition the
   * reference array.
   *
   * Calling parameters:
   *
   * reference - a K** of references sorted by its K* (x, y, z, w...) tuples
   * temporary - a K** scratch array into which to copy references
   * permutation - a vector<signed_size_t> that indications permutation of the partition coordinate
   * kdNodes - a vector<KdNode<K>>* if PREALLOCATE is defined
   *           a vector<KdNode<K>*> if PREALLOCATE is undefined
   * start - start element of the reference array
   * end - end element of the reference array
   * size - the size of the reference array
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K>* buildKdTreePresorted(K** const reference,
                                         K** const temporary,
                                         vector<signed_size_t> const& permutation,
#ifndef PREALLOCATE
                                         vector<KdNode<K>*> const& kdNodes,
#else
                                         vector<KdNode<K>>* const kdNodes,
#endif
                                         signed_size_t const start,
                                         signed_size_t const end,
                                         signed_size_t const size,
                                         signed_size_t const dim,
                                         signed_size_t const maximumSubmitDepth) {

    KdNode<K>* node = nullptr;

    // It is assumed that the reference array has been pre-sorted using the x:y:z super key.
    signed_size_t const depth = 0;

    if (end == start) {

      // Only one reference was passed to this method, so store it at this level of the tree.
      node = KdNode<K>::getKdNode(reference, kdNodes, start);
#ifdef KD_TREE_DYNAMIC_H
      node->height = 1;
#endif

} else if (end == start + 1) {

      // Two references were passed to this method in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child. 
      node = KdNode<K>::getKdNode(reference, kdNodes, start);
      node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);
#ifdef KD_TREE_DYNAMIC_H
        node->gtChild->height = 1;
        node->height = 2;
#endif

} else if (end == start + 2) {

      // Three references were passed to this method in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = KdNode<K>::getKdNode(reference, kdNodes, start + 1);
      node->ltChild = KdNode<K>::getKdNode(reference, kdNodes, start);
      node->gtChild = KdNode<K>::getKdNode(reference, kdNodes, end);
#ifdef KD_TREE_DYNAMIC_H
      node->ltChild->height = node->gtChild->height = 1;
      node->height = 2;
#endif

    } else if (end > start + 2) {

      // Four or more references were passed to this method, so use the median element of
      // the pre-sorted reference array to partition that array. Note that a < child and
      // a > child will always exist for this case.
      signed_size_t const n = end - start + 1;
      signed_size_t const median = (n + 1) >> 1;
      node = KdNode<K>::getKdNode(reference, kdNodes, median);

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so recursively build the < branch
        // of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, kdNodes, start,
                                    median - 1, size, dim, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

      } else {

        // Yes, child threads are available, so recursively build the < branch
        // of the tree with a child thread.
        auto buildFuture = async(launch::async,
                                buildKdTree,
                                reference,
                                temporary,
                                ref(permutation),
                                ref(kdNodes),
                                start,
                                median - 1,
                                size,
                                dim,
                                maximumSubmitDepth,
                                depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for build future in buildKdTreePresorted\n");
        }
      }
#ifdef KD_TREE_DYNAMIC_H
      // Compute the height at this node as the recursion unwinds.
      node->height = KdTreeDynamic<K>::computeHeight(node);
#endif

    } else if (end < start) {

      // This is an illegal condition that should never occur, so test for it last.
      ostringstream buffer;
      buffer << "\n\nerror has occurred at depth = " << depth << " : end = " << end
            << "  <  start = " << start << " in buildKdTreePresorted\n";
      throw runtime_error(buffer.str());

    }

    // Return the pointer to the root of the k-d tree.
    return node;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   * This version of createKdTree is called from the KdTreeDynamic::rebuildSubtree function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode<K>*> whose KdNodes store the (x, y, z, w...) coordinates
   * dim - the number of dimensions (required when KD_TREE_DYNAMIC_H is defined)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * p - the leading dimension
   *
   * returns: a KdTree pointer
   */
public:
  static KdTree<K>* createKdTree(vector<KdNode<K>*> const& kdNodes,
                                 size_t const dim,
                                 signed_size_t const maximumSubmitDepth,
                                 signed_size_t const p) {

    // Allocate and initialize two references arrays.
    size_t const numDimensions = dim;
    size_t const numberOfNodes = kdNodes.size();
    K*** references = new K**[2];
    for (size_t i = 0; i < 2; ++i) {
      references[i] = new K*[kdNodes.size()];
    }

    // Create a KdTree instance.
    auto tree = new KdTree<K>(numDimensions, numberOfNodes);

    // Don't allocate tuples arrays for the first references array.
    // Instead, copy pointers from the KdNode instances of the kdNodes
    // array. These pointers will be re-ordered by the KdTree::partition
    // function and then copied back to the kdNode instances by the
    // KdNode::getKdNode function. For this case where KD_TREE_DYNAMIC_H
    // is defined, the tuples will be deallocated by the ~KdNode destructor
    // when a KdNode instance is deleted from the dynamic k-d tree.
    for (size_t i = 0; i < kdNodes.size(); ++i) {
      references[0][i] = kdNodes[i]->tuple;
    }

    // For a dynamic k-d tree, it is unnecessary to sort the first references
    // array and remove duplicate coordinates, so merely specify the end index.
    signed_size_t end = kdNodes.size() - 1;

    // Determine the maximum depth of the k-d tree, which is or log2( kdNodes.size() ),
    // assuming a balanced tree.
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
    // Add the leading dimension p to the pre-computed partition coordinate (modulo
    // the number of dimensions) to permit KdTreeDynamic::balanceSubtree to build
    // a sub-tree whose root node has a non-zero partition coordinate.
    vector<signed_size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation[i] = (i + p) % numDimensions;
    }

    // Build the k-d tree with multiple threads if possible. For a dynamic k-d tree,
    // call the KdNode::buildKdTree function instead of KdNode::buildKdTreePresorted.
    tree->root = buildKdTree(references[0], references[1], permutation, kdNodes, 0, end,
                             kdNodes.size(), numDimensions, maximumSubmitDepth, 0);

    // Delete the references arrays but not the tuples arrays that they point to.
    for (size_t i = 0; i < 2; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the k-d tree.
    return tree;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   *
   * Calling parameters:
   *
   * coordinates - a vector<vector<K>> that stores the (x, y, z, w...) coordinates
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * numberOfNodes - the number of nodes counted by KdNode::verifyKdTree - returned by reference
   * allocateTime, sortTime, removeTime, kdTime, verifyTime, deallocateTime, unsortTime - execution times
   *
   * returns: a KdTree pointer
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
                                 double& deallocateTime,
                                 double& unsortTime) {

    // Allocate and initialize two references arrays.
    auto beginTime = steady_clock::now();
    size_t const numDimensions = coordinates[0].size();
    K*** references = new K**[2];
    for (size_t i = 0; i < 2; ++i) {
      references[i] = new K*[coordinates.size()];
    }

    // Create a KdTree instance.
    auto tree = new KdTree<K>(numDimensions, numberOfNodes);

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
    // Also, only the first references array is populated with tuple arrays.
    beginTime = steady_clock::now();
    MergeSort<K>::mergeSortReferenceAscending(references[0], references[1],
                                              0, coordinates.size() - 1,
                                              0, numDimensions,
                                              maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    sortTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Remove references to duplicate coordinates via one pass through the first reference array.
    beginTime = steady_clock::now();
    signed_size_t const end = KdNode<K>::removeDuplicates(references[0], 0,
                                                          numDimensions,
                                                          coordinates.size());
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

    // Start the timer to time building the k-d tree.
    beginTime = steady_clock::now();

    // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() )
    // or log2( kdNodes.size() ), depending on whether KD_TREE_DYNAMIC_H is defined,
    // and assuming a balanced tree.
    signed_size_t maxDepth = 1;
    signed_size_t size = coordinates.size();
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call of
    // the buildKdTree function because that coordinate depends only on the depth of
    // recursion, so it may be pre-computed and stored in the permutation vector.
    vector<signed_size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation[i] = i % numDimensions;
    }

    // Build the k-d tree with multiple threads if possible. For a dynamic k-d tree,
    // call the KdNode::buildKdTree function instead of KdNode::buildKdTreePresorted.

#ifndef PREALLOCATE
    tree->root = buildKdTreePresorted(references[0], references[1], permutation, kdNodes, 0, end,
                                      coordinates.size(), numDimensions, maximumSubmitDepth);
#else
    tree->root = buildKdTreePresorted(references[0], references[1], permutation, tree->kdNodes, 0, end,
                                      coordinates.size(), numDimensions, maximumSubmitDepth);
#endif

    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    kdTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Verify the k-d tree and report the number of kdNodes.
    beginTime = steady_clock::now();
    numberOfNodes = tree->root->verifyKdTree(permutation, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    verifyTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Delete the references arrays but not the tuples arrays that they point to.
    beginTime = steady_clock::now();
    for (size_t i = 0; i < 2; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // If PREALLOCATE is undefined, clear the kdNodes vector in order to measure
    // its deallocation time. Do not deallocate the KdNode instances to which it
    // points because they will be deallocated recursively by the ~KdNode destructor.

#ifndef PREALLOCATE
    kdNodes.clear();
#endif

    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    deallocateTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // The sort time is measured only for the O(kn log n) algorithm.
    sortTime = 0.0;

    // The unsort time is measured only for the Yu Cao algorithm.
    unsortTime = 0.0;

    // Return the pointer to the k-d tree.
    return tree;
  }

}; // class KdTreeNlogn

#endif // KD_TREE_NLOGN_H

