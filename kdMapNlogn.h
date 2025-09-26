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
 * 
 * -D REVERSE_NEAREST_NEIGHBORS - Enable the construction of a reverse nearest neighbors list.
 *
 * -D MEDIAN_OF_MEDIANS_CUTOFF=n - A cutoff for switching from median of medians to insertion sort
 *                                 in KdNode::partition (default 15)
 * 
 * -D MEDIAN_CUTOFF=n - A cutoff for switching from to 2 threads to calculate the median
 *                      in KdNode::partition (default 16384)
 * 
 * -D INDEX_CUTOFF=n - A cutoff for switching from to 2 threads to find the index of
 *                     the calculated median in KdNode::partition (default 512)
 * 
 * -D BIDIRECTIONAL_PARTITION - Partition an array about the median of medians proceeding
 *                              from both ends of the array instead of only the beginning.
 */

#ifndef KD_MAP_NLOGN_H
#define KD_MAP_NLOGN_H

#include "kdMapNode.h"

/* A cutoff for switching from median of medians to insertion sort in KdNode::partition */
#ifndef MEDIAN_OF_MEDIANS_CUTOFF
#define MEDIAN_OF_MEDIANS_CUTOFF 15
#endif

/* A cutoff for switching from 1 to 2 threads to calculate the median in KdNode::partition */
#ifndef MEDIAN_CUTOFF
#define MEDIAN_CUTOFF 16384
#endif

/* A cutoff for switching from 1 to 2 threads to find the index of the calculated median in KdNode::partition */
#ifndef INDEX_CUTOFF
#define INDEX_CUTOFF 1073741824 // =2^30 to disable switching to 2 threads
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
   * The swap function swaps two array elements.
   *
   * Calling parameters:
   *
   * a - KdNode** array wherein each element contains a (x,y,z,w...) tuple
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  inline
  static void swap(KdNode<K,V>** const a,
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
   * a through e - KdNode* pointers to include in the selection
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   *
   * returns a KdNode* that represents the selected KdNode
   */
private:
  inline
  static KdNode<K,V>* select_0_2(KdNode<K,V>* const a,
                                 KdNode<K,V>* const b,
                                 signed_size_t const p,
                                 signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_2(KdNode<K,V>* const a,
                                 KdNode<K,V>* const b,
                                 signed_size_t const p,
                                 signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_3_ab(KdNode<K,V>* const a,
                                    KdNode<K,V>* const b,
                                    KdNode<K,V>* const c,
                                    signed_size_t const p,
                                    signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(b->tuple, c->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_3(KdNode<K,V>* const a,
                                 KdNode<K,V>* const b,
                                 KdNode<K,V>* const c,
                                 signed_size_t const p,
                                 signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_4_ab_cd(KdNode<K,V>* const a,
                                       KdNode<K,V>* const b,
                                       KdNode<K,V>* const c,
                                       KdNode<K,V>* const d,
                                       signed_size_t const p,
                                       signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(c->tuple, a->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_4_ab(KdNode<K,V>* const a,
                                    KdNode<K,V>* const b,
                                    KdNode<K,V>* const c,
                                    KdNode<K,V>* const d,
                                    signed_size_t const p,
                                    signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(c->tuple, d->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_1_4(KdNode<K,V>* const a,
                                 KdNode<K,V>* const b,
                                 KdNode<K,V>* const c,
                                 KdNode<K,V>* const d,
                                 signed_size_t const p,
                                 signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_2_5_ab_cd(KdNode<K,V>* const a,
                                       KdNode<K,V>* const b,
                                       KdNode<K,V>* const c,
                                       KdNode<K,V>* const d,
                                       KdNode<K,V>* const e,
                                       signed_size_t const p,
                                       signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(c->tuple, a->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_2_5_ab(KdNode<K,V>* const a,
                                    KdNode<K,V>* const b,
                                    KdNode<K,V>* const c,
                                    KdNode<K,V>* const d,
                                    KdNode<K,V>* const e,
                                    signed_size_t const p,
                                    signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(c->tuple, d->tuple, p, dim) < 0) {
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
  static KdNode<K,V>* select_2_5(KdNode<K,V>* const a,
                                 KdNode<K,V>* const b,
                                 KdNode<K,V>* const c,
                                 KdNode<K,V>* const d,
                                 KdNode<K,V>* const e,
                                 signed_size_t const p,
                                 signed_size_t const dim) {
    
    if (MergeSort<K,V>::superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
   * a - a KdNode** array to recursively partition via its (x, y, z, w...) tuples array
   * start - the start index for the elements to be considered
   * n - the number of elements to consider
   * size - the size of the array of references
   * k - the element to find
   * medians - a scratch KdNode** array for the medians
   * first - the first index for the scratch array
   * p - the most significant dimension or the partition coordinate
   * dim - the number of dimensions
   * twoThreads - use two threads for median calculation
   *
   * returns - the index of the kth element in the array about which the array has been partitioned
   */
private:
  static signed_size_t partition(KdNode<K,V>** const a,
                                 signed_size_t const start,
                                 signed_size_t const n,
                                 signed_size_t const size,
                                 signed_size_t const k,
                                 KdNode<K,V>** const medians,
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
        for (j = i; j > start && MergeSort<K,V>::superKeyCompare(a[j - 1]->tuple, tmp->tuple, p, dim) > 0; --j) {
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

#ifdef DUAL_THREAD_MEDIAN
    // Is more than one thread available to calculate the medians and are
    // there sufficient medians to justify multi-threaded processing?
    //
    // NOTE, however, that NO value of MEDIAN_CUTOFF appears to improve
    // the performance of two threads relative to that of one thread.
    // Hence, the cost of spawning a child thread appears to exceed any
    // improvement in the performance of calculating the medians that
    // may be achieved via two threads.
    if (twoThreads && m > MEDIAN_CUTOFF) {

      // Yes, calculate the relative index of the middle median.
      signed_size_t mid = (m + 1) >> 1;
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
#endif
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
    // function at the expense of greater use of the superKeyCompare function.
    //
    // NOTE, however, that searching from both ends appears to degrade performance.
    signed_size_t j = n - 2;
    while (j > i) {
      if (MergeSort<K,V>::superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p, dim) < 0) {
        ++i;
      }
      else if (MergeSort<K,V>::superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p, dim) > 0) {
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
      if (MergeSort<K,V>::superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p, dim) > 0) {
        break;
      }
    }
  #else
    // Search upward from the beginning of the array a in order to minimize the use of
    // the superKeyCompare function at the expense of greater use of the swap function.
    for (signed_size_t j = 0; j < n - 1; ++j) {
      if (MergeSort<K,V>::superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p, dim) < 0) {
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
   * reference - a KdNode** array to recursively partition via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   * permutation - a vector<signed_size_t> that indications permutation of the partition coordinate
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
  static KdNode<K,V>* buildKdTree(KdNode<K,V>** reference,
                                  KdNode<K,V>** temporary,
                                  vector<signed_size_t> const& permutation,
                                  signed_size_t start,
                                  signed_size_t end,
                                  signed_size_t size,
                                  signed_size_t dim,
                                  signed_size_t maximumSubmitDepth,
                                  signed_size_t depth) {

    KdNode<K,V>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    signed_size_t const p = permutation[depth];

    if (end == start) {

      // Only one KdNode was passed to this method, so store it at this level of the tree.
      node = reference[start];
#ifdef KD_MAP_DYNAMIC_H
      node->height = 1;
#endif

    }
    else if (end == start + 1) {

      // Two KdNodes were passed to this method in unsorted order, so store the
      // start KdNode at this level of the tree and determine whether to store the
      // end KdNode as the < child or the > child.
      node = reference[start];
      if (MergeSort<K,V>::superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) > 0) {
        node->ltChild = reference[end];
#ifdef KD_MAP_DYNAMIC_H
        node->ltChild->height = 1;
        node->height = 2;
#endif
      }
      else {
        node->gtChild = reference[end];
#ifdef KD_MAP_DYNAMIC_H
        node->gtChild->height = 1;
        node->height = 2;
#endif
      }

    }
    else if (end == start + 2) {

      // Three KdNodess were passed to this method in unsorted order, so compare
      // the three KdNodes to determine which Kdnode is the median KnNode.
      // Store the median KdNode at this level of the tree, store the smallest
      // KdNode as the < child and store the largest KdNode as the > child.
      signed_size_t mid = start + 1;
      if (MergeSort<K,V>::superKeyCompare(reference[start]->tuple, reference[mid]->tuple, p, dim) < 0) {
        // reference[start] < reference[mid]
        if (MergeSort<K,V>::superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p, dim) < 0) {
          // reference[start] < reference[mid] < reference[end]
          node = reference[mid];
          node->ltChild = reference[start];
          node->gtChild = reference[end];
        }
        else {
          // reference[start] < reference[mid]; reference[end] < reference[mid]
          if (MergeSort<K,V>::superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) < 0) {
            // reference[start] < reference[end] < reference[mid]
            node = reference[end];
            node->ltChild = reference[start];
            node->gtChild = reference[mid];
          }
          else {
            // reference[end] < reference[start] < reference[mid]
            node = reference[start];
            node->ltChild = reference[end];
            node->gtChild = reference[mid];
          }
        }
      }
      else {
        // reference[mid] < reference[start]
        if (MergeSort<K,V>::superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) < 0) {
          // reference[mid] < reference[start] < reference[end]
          node = reference[start];
          node->ltChild = reference[mid];
          node->gtChild = reference[end];
        }
        else {
          // reference[mid] < reference[start]; reference[end] < reference[start]
          if (MergeSort<K,V>::superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p, dim) < 0) {
            // reference[mid] < reference[end] < reference[start]
            node = reference[end];
            node->ltChild = reference[mid];
            node->gtChild = reference[start];
          }
          else {
            // reference[end] < reference[mid] < reference[start]
            node = reference[mid];
            node->ltChild = reference[end];
            node->gtChild = reference[start];
          }
        }
      }

#ifdef KD_MAP_DYNAMIC_H
      node->ltChild->height = node->gtChild->height = 1;
      node->height = 2;
#endif

    }
    else if (end > start + 2) {

      // Four or more references were passed to this method, so calculate the offset
      // of the median element. partition the reference array about its median element,
      // which is the kth element as calculated below.
      signed_size_t const n = end - start + 1;
      signed_size_t const k = (n + 1) >> 1;

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so find the median element then
        // partition the reference array about it.  Store the median element
        // from the reference array in a new KdNode.
        signed_size_t const median = partition(reference, start, n, size, k, temporary, start, p, dim, false);
        node = reference[median];

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, start,
                                    median - 1, size, dim, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so find the median element then partition
        // the reference array about it.  Store the median element from the reference
        // array in a new KdNode.
        signed_size_t median = partition(reference, start, n, size, k, temporary, start, p, dim, true);
        node = reference[median];

        // Recursively build the < branch of the tree with a child thread.
        auto buildFuture = async(launch::async, buildKdTree,
                                reference,
                                temporary,
                                ref(permutation),
                                start,
                                median - 1,
                                size,
                                dim,
                                maximumSubmitDepth,
                                depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

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
 * The buildKdTreePresorted function method builds a k-d tree by using
 * the median of the pre-sorted reference array to partition that array,
 * then calls the buildKdTree function to recursively partition the
 * reference array.
 *
 * Calling parameters:
 *
 * reference - a KdNode** array sorted via its (x, y, z, w...) tuples array
 * temporary - a KdNode** temporary array from which to copy sorted results;
 * permutation - a vector<signed_size_t> that indications permutation of the partition coordinate
 * start - start element of the reference array
 * end - end element of the reference array
 * size - the size of the reference array
 * dim - the number of dimensions
 * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
 *
 * returns a KdNode pointer to the root of the k-d tree
 */
private:
  static KdNode<K,V>* buildKdTreePresorted(KdNode<K,V>** reference,
                                          KdNode<K,V>** temporary,
                                          vector<signed_size_t> const& permutation,
                                          signed_size_t start,
                                          signed_size_t end,
                                          signed_size_t size,
                                          signed_size_t dim,
                                          signed_size_t maximumSubmitDepth) {

    KdNode<K,V>* node = nullptr;

    // It is assumed that the reference array has been pre-sorted using the x:y:z:w... super key.
    signed_size_t const depth = 0;

    if (end == start) {

      // Only one KdNode was passed to this method, so store it at this level of the tree.
      node = reference[start];
#ifdef KD_MAP_DYNAMIC_H
      node->height = 1;
#endif

    }
    else if (end == start + 1) {

      // Two KdNodes were passed to this method in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child. 
      node = reference[start];
      node->gtChild = reference[end];
#ifdef KD_MAP_DYNAMIC_H
        node->gtChild->height = 1;
        node->height = 2;
#endif

    }
    else if (end == start + 2) {

      // Three KdNodes were passed to this method in sorted order, so
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

      // Four or more KdNodes were passed to this method, so use the median element of
      // the pre-sorted reference array to partition the reference array.
      signed_size_t const n = end - start + 1;
      signed_size_t const median = (n + 1) >> 1;
      node = reference[median];

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so recursively build the < branch
        // of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, start, median - 1,
                                    size, dim, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1, end,
                                    size, dim, maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so recursively build the < branch
        // of the tree with a child thread. The recursive call to buildKdTree
        // must be placed in a lambda expression because buildKdTree is a template
        // not a function.
        auto buildFuture = async(launch::async, buildKdTree,
                                reference,
                                temporary,
                                ref(permutation),
                                start,
                                median - 1,
                                size,
                                dim,
                                maximumSubmitDepth,
                                depth + 1);

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for build future in buildKdTreePresorted\n");
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
   * dim - the number of dimensions (required when KD_MAP_DYNAMIC_H is defined)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
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

    // Allocate two references arrays.
    size_t numDimensions = dim;
    KdNode<K,V>*** references = new KdNode<K,V>**[2];
    for (size_t i = 0; i < 2; ++i) {
      references[i] = new KdNode<K,V>*[kdNodes.size()];
    }

    // Don't allocate KdNode instances for the first references array.
    // Instead, copy pointers from the KdNode instances of the kdNodes
    // vector. These pointers will be re-ordered by the KdTree::partition
    // function. For this case where KD_MAP_DYNAMIC_H is defined, the
    // tuples will be deallocated by the ~KdNode destructor when a
    // KdNode instance is deleted from the dynamic k-d tree.
    for (size_t i = 0; i < kdNodes.size(); ++i) {
      references[0][i] = kdNodes[i];
    }

    // For a dynamic k-d tree, it is unnecessary to sort the first references
    // array and remove duplicate coordinates, so merely specify the end index.
    signed_size_t end = kdNodes.size() - 1;

    // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() )
    // or log2( kdNodes.size() ), depending on whether KD_MAP_DYNAMIC_H is defined,
    // and assuming a balanced tree.
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
    tree->root = buildKdTree(references[0], references[1], permutation, 0, end,
                             kdNodes.size(), numDimensions, maximumSubmitDepth, 0);

    // Delete the references arrays but not the KdNodes instances that they point to
    // because those KdNodes instances will be deleted by the ~KdTree destructor.
    for (size_t i = 0; i < 2; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the KdTree instance.
    return tree;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTreePresorted function.
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

    // Allocate two references arrays.
    size_t numDimensions = coordinates[0].first.size();

    auto beginTime = steady_clock::now();
    KdNode<K,V>*** references = new KdNode<K,V>**[2];
    for (size_t i = 0; i < 2; ++i) {
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
    MergeSort<K,V>::mergeSortReferenceAscending(references[0], references[1],
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

    // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() )
    // or log2( kdNodes.size() ), depending on whether KD_MAP_DYNAMIC_H is defined,
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

    // Build the k-d tree with multiple threads if possible.
    tree->root = buildKdTreePresorted(references[0], references[1], permutation, 0, end,
                                      coordinates.size(), numDimensions, maximumSubmitDepth);
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
    numberOfNodes = tree->root->verifyKdTree(permutation, numDimensions, maximumSubmitDepth, 0);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    verifyTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
  
    // Delete the references arrays but not the KdNodes instances that they point to
    // because those KdNodes instances will be deleted by the ~KdTree destructor.
    beginTime = steady_clock::now();
    for (size_t i = 0; i < 2; ++i) {
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

#endif // KD_MAP_NLOGN_H
