/*
 * Copyright (c) 2015, 2021 Russell A. Brown
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
 * The k-d tree was described by Jon Bentley in "Multidimensional Binary Search Trees
 * Used for Associative Searching", CACM 18(9): 509-517, 1975.  For k dimensions and
 * n elements of data, a balanced k-d tree is built in O(n log n) time by finding the
 * median of the data at each level of the tree via the "median of medians" algorithm
 * described by Manuel Blum, et al. in "Time Bounds for Selection", Journal of Computer
 * and System Sciences, 7: 448-461, 1973.
 *
 * Gnu g++ compilation options are: -lm -O3 -DTEST_KD_TREE
 *
 * Optional compilation options are:
 *
 * -D K_DIMENSIONALITY=n - The k-dimensionality of the k-d tree (default 3).
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * -D MEDIAN_OF_MEDIANS_CUTOFF=n - A cutoff for switching from median of medians to insertion sort
 *                                 in KdNode::partition (default 15)
 * -D MEDIAN_CUTOFF=n - A cutoff for switching from to 2 threads to calculate the median
                        in KdNode::partition (default 16384)
 * -D INDEX_CUTOFF=n - A cutoff for switching from to 2 threads to find the index of
                       the calculated median in KdNode::partition (default 512)
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 * -D DUAL_THREAD_MEDIAN - Calculate the medians with two threads.
 * -D DUAL_THREAD_INDEX - Find the index of the median of medians with two threads.
 * -D BIDIRECTIONAL_PARTITION - Partition an array about the median of medians proceeding
 *                              from both ends of the array instead of only the beginning.
 * -D MACH - Use a Mach equivalent to the clock_gettime(CLOCK_REALTIME, &time) function
 *           but this option appears to no longer be necessary.
 */

#include <limits.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <exception>
#include <future>

/* Here is the default k-dimensionality of the k-d tree. */
#ifndef K_DIMENSIONALITY
#define K_DIMENSIONALITY 3
#endif

/* A cutoff for switching from merge sort to insertion sort in the KdNode::mergeSort* functions */
#ifndef INSERTION_SORT_CUTOFF
#define INSERTION_SORT_CUTOFF 15
#endif

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
#define INDEX_CUTOFF = 512
#endif

/*
 * This type is the signed equivalent of size_t and might be equivalent to intmax_t
 */
typedef std::streamsize signed_size_t;

/*
 * These are the types used for the test. Change the intrisic types in
 * these typedefs to test the k-d tree with different intrisic types.
 */
typedef int64_t kdKey_t;
typedef signed_size_t kdValue_t;

/*
 * Create an alternate to clock_gettime(CLOCK_REALTIME, &time) for Mach. See
 * http://stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x
 * However, it appears that later versions of Mac OS X support clock_gettime(),
 * so this alternative may no longer be necessary for Mac OS X.
 */
#ifdef MACH
#include <mach/mach_time.h>

#define MACH_NANO (+1.0E-9)
#define MACH_GIGA UINT64_C(1000000000)

static double mach_timebase = 0.0;
static uint64_t mach_timestart = 0;

struct timespec getTime(void) {
  // be more careful in a multithreaded environement
  if (!mach_timestart) {
    mach_timebase_info_data_t tb = { 0, 1 }; // Initialize tb.numer and tb.denom
    mach_timebase_info(&tb);
    mach_timebase = tb.numer;
    mach_timebase /= tb.denom;
    mach_timestart = mach_absolute_time();
  }
  struct timespec t;
  double diff = (mach_absolute_time() - mach_timestart) * mach_timebase;
  t.tv_sec = diff * MACH_NANO;
  t.tv_nsec = diff - (t.tv_sec * MACH_GIGA);
  return t;
}
#else

#if defined(_WIN32) || defined(_WIN64)
//see https://stackoverflow.com/questions/5404277/porting-clock-gettime-to-windows/5404467#5404467
#define NOMINMAX // Prevent Windows from getting confused about std::min vs. min, etc.
#include <windows.h>

int clock_gettime(int, struct timespec* tv)
{
  static int initialized = 0;
  static LARGE_INTEGER freq, startCount;
  static struct timespec tv_start;
  LARGE_INTEGER curCount;
  time_t sec_part;
  long nsec_part;

  if (!initialized) {
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&startCount);
    timespec_get(&tv_start, TIME_UTC);
    initialized = 1;
  }

  QueryPerformanceCounter(&curCount);

  curCount.QuadPart -= startCount.QuadPart;
  sec_part = curCount.QuadPart / freq.QuadPart;
  nsec_part = (long)((curCount.QuadPart - (sec_part * freq.QuadPart))
                     * 1000000000UL / freq.QuadPart);

  tv->tv_sec = tv_start.tv_sec + sec_part;
  tv->tv_nsec = tv_start.tv_nsec + nsec_part;
  if (tv->tv_nsec >= 1000000000UL) {
    tv->tv_sec += 1;
    tv->tv_nsec -= 1000000000UL;
  }
  return 0;
}

#define CLOCK_REALTIME 0

#endif

struct timespec getTime(void) {
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  return time;
}
#endif

/* A forward reference to NearestNeighborHeap */
template <typename, typename, signed_size_t>
class NearestNeighborHeap;

/* One node of a k-d tree where K is key type, V is value type, and N is tuple dimensionality */
template <typename K, typename V, signed_size_t N>
class KdNode {
public:
  K tuple[N];
private:
  KdNode<K,V,N>* ltChild;
  KdNode<K,V,N>* gtChild;
  KdNode<K,V,N>* duplicates;
  V value;

public:
  KdNode(V const value) { // Pass non-primitive types as 'V const&'
    this->ltChild = this->gtChild = this->duplicates = nullptr; // redundant
    this->value = value;
  }

public:
  ~KdNode() {
    // Delete each KdNode from the duplicates list.
    auto nextPtr = this->duplicates;
    while (nextPtr != nullptr) {
      auto tempPtr = nextPtr;
      nextPtr = nextPtr->duplicates;
      delete tempPtr;
    }
  }

public:
  K const* getTuple() {
    return this->tuple;
  }

  /*
   * The superKeyCompare function compares two K arrays in all k dimensions,
   * and uses the sorting or partition coordinate as the most significant dimension.
   *
   * Calling parameters:
   *
   * a - a K*
   * b - a K*
   * p - the most significant dimension
   *
   * returns a K result of comparing two K arrays
   */
private:
  inline
  static K superKeyCompare(K const* a, K const* b, signed_size_t p) {
    // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
    K diff = a[p] - b[p];
    for (signed_size_t i = 1; diff == 0 && i < N; i++) {
      signed_size_t r = i + p;
      // A fast alternative to the modulus operator for (i + p) < 2 * N.
      r = (r < N) ? r : r - N;
      diff = a[r] - b[r];
    }
    return diff;
  }

  /*
   * The following four merge sort functions are adapted from the mergesort function that is shown
   * on p. 166 of Robert Sedgewick's "Algorithms in C++", Addison-Wesley, Reading, MA, 1992.
   * That elegant implementation of the merge sort algorithm eliminates the requirement to test
   * whether the upper and lower halves of an auxiliary array have become exhausted during the
   * merge operation that copies from the auxiliary array to a result array.  This elimination is
   * made possible by inverting the order of the upper half of the auxiliary array and by accessing
   * elements of the upper half of the auxiliary array from highest address to lowest address while
   * accessing elements of the lower half of the auxiliary array from lowest address to highest
   * address.
   *
   * The following four merge sort functions also implement two suggestions from p. 275 of Robert
   * Sedgewick's and Kevin Wayne's "Algorithms 4th Edition", Addison-Wesley, New York, 2011.  The
   * first suggestion is to replace merge sort with insertion sort when the size of the array to
   * sort falls below a threshold.  The second suggestion is to avoid unnecessary copying to the
   * auxiliary array prior to the merge step of the algorithm by implementing two versions of
   * merge sort and by applying some "recursive trickery" to arrange that the required result is
   * returned in an auxiliary array by one version and in a result array by the other version.
   * The following four merge sort methods build upon this suggestion and return their result in
   * either ascending or descending order, as discussed on pp. 173-174 of Robert Sedgewick's
   * "Algorithms in C++", Addison-Wesley, Reading, MA, 1992.
   *
   * During multi-threaded execution, the upper and lower halves of the result array may be filled
   * from the auxiliary array (or vice versa) simultaneously by two threads.  The lower half of the
   * result array is filled by accessing elements of the upper half of the auxiliary array from highest
   * address to lowest address while accessing elements of the lower half of the auxiliary array from
   * lowest address to highest address, as explained above for elimination of the test for exhaustion.
   * The upper half of the result array is filled by addressing elements from the upper half of the
   * auxiliary array from lowest address to highest address while accessing the elements from the lower
   * half of the auxiliary array from highest address to lowest address.  Note: for the upper half
   * of the result array, there is no requirement to test for exhaustion provided that the upper half
   * of the result array never comprises more elements than the lower half of the result array.  This
   * provision is satisfied by computing the median address of the result array as shown below for
   * all four merge sort methods.
   *
   *
   * The mergeSortReferenceAscending function recursively subdivides the reference array then
   * merges the elements in ascending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   n  * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceAscending(KdNode<K,V,N>** reference, KdNode<K,V,N>** temporary,
                                          signed_size_t low, signed_size_t high, signed_size_t p,
                                          signed_size_t maximumSubmitDepth, signed_size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in ascending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) < 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in ascending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortTemporaryAscending, reference, temporary,
                                                  low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in ascending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
                                           for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                             reference[k] =
                                               (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) <= 0)
                                               ? temporary[i++] : temporary[j--];
                                           }
                                         });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in ascending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) > 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        KdNode<K,V,N>* tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p) > 0; --j) {
          reference[j] = reference[j - 1];
        }
        reference[j] = tmp;
      }
    }
  }

  /*
   * The mergeSortReferenceDecending function recursively subdivides the reference array then
   * merges the elements in descending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceDescending(KdNode<K,V,N>** reference, KdNode<K,V,N>** temporary,
                                           signed_size_t low, signed_size_t high, signed_size_t p,
                                           signed_size_t maximumSubmitDepth, signed_size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in descending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) > 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in descending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortTemporaryDescending, reference, temporary,
                                                  low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in descending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
                                           for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                             reference[k] =
                                               (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) >= 0)
                                               ? temporary[i++] : temporary[j--];
                                           }
                                         });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in descending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) < 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        KdNode<K,V,N>* tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p) < 0; --j) {
          reference[j] = reference[j - 1];
        }
        reference[j] = tmp;
      }
    }
  }

  /*
   * The mergeSortTemporaryAscending function recursively subdivides the reference array then
   * merges the elements in ascending order and leaves the result in the temporary array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryAscending(KdNode<K,V,N>** reference, KdNode<K,V,N>** temporary,
                                          signed_size_t low, signed_size_t high, signed_size_t p,
                                          signed_size_t maximumSubmitDepth, signed_size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in ascending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) < 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in ascending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortReferenceAscending, reference, temporary,
                                                  low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in ascending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
                                           for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                             temporary[k] =
                                               (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) <= 0)
                                               ? reference[i++] : reference[j--];
                                           }
                                         });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in ascending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) > 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is John Robinson's implementation of insertion sort that sorts in ascending order
      // and leaves the result in the temporary array.
      temporary[high] = reference[high];
      signed_size_t i;
      signed_size_t j; // MUST be signed because it can decrement to -1
      for (j = high - 1; j >= low; --j) {
        for (i = j; i < high; ++i) {
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p) > 0) {
            temporary[i] = temporary[i + 1];
          }
          else {
            break;
          }
        }
        temporary[i] = reference[j];
      }
    }
  }

  /*
   * The mergeSortTemporaryDecending function recursively subdivides the reference array
   * then merges the elements in descending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryDescending(KdNode<K,V,N>** reference, KdNode<K,V,N>** temporary,
                                           signed_size_t low, signed_size_t high, signed_size_t p,
                                           signed_size_t maximumSubmitDepth, signed_size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in descending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) > 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in descending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortReferenceDescending, reference, temporary,
                                                  low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in descending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
                                           for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                             temporary[k] =
                                               (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) >= 0)
                                               ? reference[i++] : reference[j--];
                                           }
                                         });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in descending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) < 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is John Robinson's implementation of insertion sort that sorts in descending order
      // and leaves the result in the temporary array.
      temporary[high] = reference[high];
      signed_size_t i;
      signed_size_t j; // MUST be signed because it can decrement to -1
      for (j = high - 1; j >= low; --j) {
        for (i = j; i < high; ++i) {
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p) < 0) {
            temporary[i] = temporary[i + 1];
          }
          else {
            break;
          }
        }
        temporary[i] = reference[j];
      }
    }
  }

  /*
   * The removeDuplicates function checks the validity of the merge sort and
   * removes duplicates from the kdNodes array.
   *
   * Calling parameters:
   *
   * kdNodes - a KdNode** array that has been sorted via merge sort according to (x,y,z,w...) tuples
   * i - the leading dimension for the super key
   *
   * returns the end index of the reference array following removal of duplicate elements
   */
private:
  inline
  static signed_size_t removeDuplicates(KdNode<K,V,N>** kdNodes, signed_size_t i, signed_size_t size) {
    signed_size_t end = 0;
    for (signed_size_t j = 1; j < size; ++j) {
      K compare = superKeyCompare(kdNodes[j]->tuple, kdNodes[end]->tuple, i);
      if (compare < 0) {
        std::cout << "merge sort failure: superKeyCompare(kdNodes[" << j << "], kdNodes["
                  << j - 1 << "], (" << i << ") = " << compare << std::endl;
        exit(1);
      }
      else if (compare > 0) {
        // Keep the jth element of the kdNodes array.
        kdNodes[++end] = kdNodes[j];
      } else {
        // Discard the jth element of the kdNodes array and prepend it to the duplicates list.
        kdNodes[j]->duplicates = kdNodes[end]->duplicates;
        kdNodes[end]->duplicates = kdNodes[j];
      }
    }
    return end;
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
  static void swap(KdNode<K,V,N>** a, signed_size_t i, signed_size_t j) {
    KdNode<K,V,N>* t = a[i];
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
   *
   * returns a KdNode* that represents the selected KdNode
   */
private:
  inline
  static KdNode<K,V,N>* select_0_2(KdNode<K,V,N>* a,
                                   KdNode<K,V,N>* b,
                                   signed_size_t p) {
    if (superKeyCompare(a->tuple, b->tuple, p) < 0) {
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
  static KdNode<K,V,N>* select_1_2(KdNode<K,V,N>* a,
                                   KdNode<K,V,N>* b,
                                   signed_size_t p) {
    if (superKeyCompare(a->tuple, b->tuple, p) < 0) {
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
  static KdNode<K,V,N>* select_1_3_ab(KdNode<K,V,N>* a,
                                      KdNode<K,V,N>* b,
                                      KdNode<K,V,N>* c,
                                      signed_size_t p) {
    if (superKeyCompare(b->tuple, c->tuple, p) < 0) {
      // a < b < c
      return b;
    }
    else {
      // a ? c < b
      return select_1_2(a, c, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_1_3(KdNode<K,V,N>* a,
                                   KdNode<K,V,N>* b,
                                   KdNode<K,V,N>* c,
                                   signed_size_t p) {
    if (superKeyCompare(a->tuple, b->tuple, p) < 0) {
      // a < b
      return select_1_3_ab(a, b, c, p);
    }
    else {
      // b < a
      return select_1_3_ab(b, a, c, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_1_4_ab_cd(KdNode<K,V,N>* a,
                                         KdNode<K,V,N>* b,
                                         KdNode<K,V,N>* c,
                                         KdNode<K,V,N>* d,
                                         signed_size_t p) {
    if (superKeyCompare(c->tuple, a->tuple, p) < 0) {
      // c < a < b && a ? d so c is eliminated and a ? d
      return select_0_2(a, d, p);
    }
    else {
      // a < b ? c < d so a is eliminated and b ? c
      return select_0_2(b, c, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_1_4_ab(KdNode<K,V,N>* a,
                                      KdNode<K,V,N>* b,
                                      KdNode<K,V,N>* c,
                                      KdNode<K,V,N>* d,
                                      signed_size_t p) {
    if (superKeyCompare(c->tuple, d->tuple, p) < 0) {
      // a < b && c < d
      return select_1_4_ab_cd(a, b, c, d, p);
    }
    else {
      // a < b && d < c
      return select_1_4_ab_cd(a, b, d, c, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_1_4(KdNode<K,V,N>* a,
                                   KdNode<K,V,N>* b,
                                   KdNode<K,V,N>* c,
                                   KdNode<K,V,N>* d,
                                   signed_size_t p) {
    if (superKeyCompare(a->tuple, b->tuple, p) < 0) {
      // a < b
      return select_1_4_ab(a, b, c, d, p);
    }
    else {
      // b < a
      return select_1_4_ab(b, a, c, d, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_2_5_ab_cd(KdNode<K,V,N>* a,
                                         KdNode<K,V,N>* b,
                                         KdNode<K,V,N>* c,
                                         KdNode<K,V,N>* d,
                                         KdNode<K,V,N>* e,
                                         signed_size_t p) {
    if (superKeyCompare(c->tuple, a->tuple, p) < 0) {
      // c < a < b && c < d ? e so c is eliminated and a < b && d ? e
      return select_1_4_ab(a, b, d, e, p);
    }
    else {
      // a < b ? c && c < d ? e && b ? e so a is eliminated and c < d && b ? e
      return select_1_4_ab(c, d, b, e, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_2_5_ab(KdNode<K,V,N>* a,
                                      KdNode<K,V,N>* b,
                                      KdNode<K,V,N>* c,
                                      KdNode<K,V,N>* d,
                                      KdNode<K,V,N>* e,
                                      signed_size_t p) {
    if (superKeyCompare(c->tuple, d->tuple, p) < 0) {
      // a < b && c < d
      return select_2_5_ab_cd(a, b, c, d, e, p);
    }
    else {
      // a < b && d < c
      return select_2_5_ab_cd(a, b, d, c, e, p);
    }
  }

private:
  inline
  static KdNode<K,V,N>* select_2_5(KdNode<K,V,N>* a,
                                   KdNode<K,V,N>* b,
                                   KdNode<K,V,N>* c,
                                   KdNode<K,V,N>* d,
                                   KdNode<K,V,N>* e,
                                   signed_size_t p) {
    if (superKeyCompare(a->tuple, b->tuple, p) < 0) {
      // a < b
      return select_2_5_ab(a, b, c, d, e, p);
    }
    else {
      // b < a
      return select_2_5_ab(b, a, c, d, e, p);
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
   * twoThreads - use two threads for median calculation
   *
   * returns - the index of the kth element in the array about which the array has been partitioned
   */
private:
  static signed_size_t partition(KdNode<K,V,N>** a, signed_size_t start, signed_size_t n, signed_size_t size,
                                 signed_size_t k, KdNode<K,V,N>** medians, signed_size_t first, signed_size_t p,
                                 bool twoThreads) {

    if (n <= 0 || n > size) {
      std::cout << "Error in n = " << n << "  size = " << size << std::endl;
    }
    if (k <= 0 || k > n) {
      std::cout << "Error in k = " << k << std::endl;
    }
    if (start + n > size) {
      std::cout << "Error in start = " << start << "  n = " << n << "  size = " << size << std::endl;
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
        KdNode<K,V,N>* tmp = a[i];
        signed_size_t j;
        for (j = i; j > start && superKeyCompare(a[j - 1]->tuple, tmp->tuple, p) > 0; --j) {
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
      std::future<void> medianFuture =
        std::async(std::launch::async, [&] {
                                         for (signed_size_t firstOfGroup = 0, i = 0; i < mid; ++i) {

                                           // Find the median of the group of GROUP_SIZE elements via select_2_5.
                                           medians[first + i] = select_2_5(a[start + firstOfGroup],
                                                                           a[start + firstOfGroup + 1],
                                                                           a[start + firstOfGroup + 2],
                                                                           a[start + firstOfGroup + 3],
                                                                           a[start + firstOfGroup + 4],
                                                                           p);

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
                                        p);

        // Update the index to the next group of GROUP_SIZE elements.
        startOfGroup += GROUP_SIZE;
      }

      // Wait for the child thread to finish execution.
      try {
        medianFuture.get();
      }
      catch (std::exception const& e) {
        std::cout << "caught exception " << e.what() << std::endl;
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
                                        p);

        // Update the index to the next group of GROUP_SIZE elements.
        startOfGroup += GROUP_SIZE;
      }
    }

    // Calculate and check the number of remaining elements.
    signed_size_t remainingElements = n - startOfGroup;
    if (remainingElements < 0 || remainingElements >= GROUP_SIZE) {
      std::cout << "Error: incorrect group calculation";
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
                                        p);
        ++m;
        break;
      case 3:
        medians[first + m] = select_1_3(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        p);
        ++m;
        break;
      case 4:
        medians[first + m] = select_1_4(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        p);
        ++m;
        break;
      default:
        std::cout << "Error: unhandled case in switch: remainingElements = " << remainingElements << std::endl;
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
    KdNode<K,V,N> const* medianOfMedians =
      medians[partition(medians, first, m, first + m, (m + 1) >> 1, medians, first + m, p, twoThreads)];

    // Find the index of the median of medians and swap it into a[start + n - 1]
    // so that it is not examined during partitioning of the array a.
    //
    // Is more than one thread available to find the index of the median of medians
    // and are there sufficient array elements to justify dual-threaded processing?
    //
    // NOTE, however, that NO value of INDEX_CUTOFF appears to improve
    // the performance of two threads relative to that of one thread.
    // Hence, the cost of spawning a child thread appears to exceed any
    // improvement in the performance of finding the index of the median
    // of medians that may be achieved via two threads.
#ifdef DUAL_THREAD_INDEX
    if (twoThreads && n > INDEX_CUTOFF) {

      // Yes, more than one thread is available, so calculate the relative index of the middle element.
      signed_size_t middle = (n + 1) >> 1;

      // Search for the index in the lower half of the array a with a child thread.
      std::future<void> indexFuture =
        std::async(std::launch::async, [&] {
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
      catch (std::exception const& e) {
        std::cout << "caught exception " << e.what() << std::endl;
      }
    }
    else
#endif
    {
      // No, only one thread is available to find the index of the median of medians.
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
      if (superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p) < 0) {
        ++i;
      }
      else if (superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p) > 0) {
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
      if (superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p) > 0) {
        break;
      }
    }
#else
    // Search upward from the beginning of the array a in order to minimize the use of
    // the superKeyCompare function at the expense of greater use of the swap function.
    for (signed_size_t j = 0; j < n - 1; ++j) {
      if (superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p) < 0) {
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
      return partition(a, start, i, size, k, medians, first, p, twoThreads);

    }
    else if (k > i + 1) {

      // The median of medians occupies a position above i, so partition
      // the array elements of the > subset; for this subset, the
      // original kth element is not the kth element of this subset
      // because i + 1 elements are in the < subset.
      return partition(a, start + i + 1, n - i - 1, size, k - i - 1,
                       medians, first, p, twoThreads);

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
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K,V,N>* buildKdTree(KdNode<K,V,N>** reference,
                                    KdNode<K,V,N>** temporary,
                                    std::vector<signed_size_t> const& permutation,
                                    signed_size_t start, signed_size_t end, signed_size_t size,
                                    signed_size_t maximumSubmitDepth,
                                    signed_size_t depth) {

    KdNode<K,V,N>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    signed_size_t p = permutation[depth];

    if (end == start) {

      // Only one KdNode was passed to this method, so store it at this level of the tree.
      node = reference[start];

    }
    else if (end == start + 1) {

      // Two KdNodes were passed to this method in unsorted order, so store the
      // start KdNode at this level of the tree and determine whether to store the
      // end KdNode as the < child or the > child.
      node = reference[start];
      if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p) > 0) {
        node->ltChild = reference[end];
      }
      else {
        node->gtChild = reference[end];
      }

    }
    else if (end == start + 2) {

      // Three KdNodess were passed to this method in unsorted order, so compare
      // the three KdNodes to determine which Kdnode is the median KnNode.
      // Store the median KdNode at this level of the tree, store the smallest
      // KdNode as the < child and store the largest KdNode as the > child.
      signed_size_t mid = start + 1;
      if (superKeyCompare(reference[start]->tuple, reference[mid]->tuple, p) < 0) {
        // reference[start] < reference[mid]
        if (superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p) < 0) {
          // reference[start] < reference[mid] < reference[end]
          node = reference[mid];
          node->ltChild = reference[start];
          node->gtChild = reference[end];
        }
        else {
          // reference[start] < reference[mid]; reference[end] < reference[mid]
          if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p) < 0) {
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
        if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p) < 0) {
          // reference[mid] < reference[start] < reference[end]
          node = reference[start];
          node->ltChild = reference[mid];
          node->gtChild = reference[end];
        }
        else {
          // reference[mid] < reference[start]; reference[end] < reference[start]
          if (superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p) < 0) {
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

    }
    else if (end > start + 2) {

      // Four or more references were passed to this method, so calculate the offset
      // of the median element. partition the reference array about its median element,
      // which is the kth element as calculated below.
      signed_size_t n = end - start + 1;
      signed_size_t k = (n + 1) >> 1;

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so find the median element then
        // partition the reference array about it.  Store the median element
        // from the reference array in a new KdNode.
        signed_size_t median = partition(reference, start, n, size, k, temporary, start, p, false);
        node = reference[median];

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, start,
                                    median - 1, size, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so find the median element then partition
        // the reference array about it.  Store the median element from the reference
        // array in a new KdNode.
        signed_size_t median = partition(reference, start, n, size, k, temporary, start, p, true);
        node = reference[median];

        // Recursively build the < branch of the tree with a child thread.
        // The recursive call to buildKdTree must be placed in a lambda
        // expression because buildKdTree is a template not a function.
        std::future<KdNode<K,V,N>*> buildFuture =
          std::async(std::launch::async, [&] {
                                           return buildKdTree(reference, temporary, permutation, start,
                                                              median - 1, size, maximumSubmitDepth, depth + 1);
                                         });

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else if (end < start) {

      // This is an illegal condition that should never occur, so test for it last.
      std::cout << "error has occurred at depth = " << depth << " : end = " << end
                << "  <  start = " << start << std::endl;
      exit(1);

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
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K,V,N>* buildKdTreePresorted(KdNode<K,V,N>** reference,
                                             KdNode<K,V,N>** temporary,
                                             std::vector<signed_size_t> const& permutation,
                                             signed_size_t start, signed_size_t end, signed_size_t size,
                                             signed_size_t maximumSubmitDepth) {

    KdNode<K,V,N>* node = nullptr;

    // It is assumed that the reference array has been pre-sorted using the x:y:z:w... super key.
    signed_size_t depth = 0;

    if (end == start) {

      // Only one KdNode was passed to this method, so store it at this level of the tree.
      node = reference[start];

    }
    else if (end == start + 1) {

      // Two KdNodes were passed to this method in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child. 
      node = reference[start];
      node->gtChild = reference[end];

    }
    else if (end == start + 2) {

      // Three KdNodes were passed to this method in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = reference[start + 1];
      node->ltChild = reference[start];
      node->gtChild = reference[end];

    }
    else if (end > start + 2) {

      // Four or more KdNodes were passed to this method, so use the median element of
      // the pre-sorted reference array to partition the reference array.
      signed_size_t n = end - start + 1;
      signed_size_t median = (n + 1) >> 1;
      node = reference[median];

      // Build the < branch of the tree with a child thread at as many levels of the
      // tree as possible.  Create the child thread as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so recursively build the < branch
        // of the tree with the current thread.
        node->ltChild = buildKdTree(reference, temporary, permutation, start, median - 1,
                                    size, maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1, end,
                                    size, maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so recursively build the < branch
        // of the tree with a child thread. The recursive call to buildKdTree
        // must be placed in a lambda expression because buildKdTree is a template
        // not a function.
        std::future<KdNode<K,V,N>*> buildFuture =
          std::async(std::launch::async, [&] {
                                           return buildKdTree(reference, temporary, permutation, start,
                                                              median - 1, size, maximumSubmitDepth, depth + 1);
                                         });

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(reference, temporary, permutation, median + 1,
                                    end, size, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else if (end < start) {

      // This is an illegal condition that should never occur, so test for it last.
      std::cout << "error has occurred at depth = " << depth << " : end = " << end
                << "  <  start = " << start << std::endl;
      exit(1);

    }

    // Return the pointer to the root of the k-d tree.
    return node;
  }

  /*
   * The verifyKdTree function walks the k-d tree and checks that the
   * children of a node are in the correct branch of that node.
   *
   * Calling parameters:
   *
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   *
   * returns: a count of the number of kdNodes in the k-d tree
   */
private:
  signed_size_t verifyKdTree(signed_size_t maximumSubmitDepth, signed_size_t depth) {

    signed_size_t count = 1;

    // The partition cycles as x, y, z, w...
    signed_size_t p = depth % N;

    if (ltChild != nullptr) {
      if (ltChild->tuple[p] > tuple[p]) {
        std::cout << "child is > node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(ltChild->tuple, tuple, p) >= 0) {
        std::cout << "child is >= node!" << std::endl;
        exit(1);
      }
    }
    if (gtChild != nullptr) {
      if (gtChild->tuple[p] < tuple[p]) {
        std::cout << "child is < node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(gtChild->tuple, tuple, p) <= 0) {
        std::cout << "child is <= node" << std::endl;
        exit(1);
      }
    }

    // Verify the < branch with a child thread at as many levels of the tree as possible.
    // Create the child thread as high in the tree as possible for greater utilization.

    // Is a child thread available to build the < branch?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so verify the < branch with the current thread.
      if (ltChild != nullptr) {
        count += ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Then verify the > branch with the current thread.
      if (gtChild != nullptr) {
        count += gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }
    }
    else {

      // Yes, so verify the < branch with a child thread. Note that a
      // lambda is required to instantiate the verifyKdTree template.
      std::future<signed_size_t> verifyFuture;
      if (ltChild != nullptr) {
        verifyFuture =
          std::async(std::launch::async, [&] {
                                           return ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
                                         });
      }

      // And simultaneously verify the > branch with the current thread.
      signed_size_t gtCount = 0;
      if (gtChild != nullptr) {
        gtCount = gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Wait for the child thread to finish execution.
      signed_size_t ltCount = 0;
      if (ltChild != nullptr) {
        try {
          ltCount = verifyFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }
      count += ltCount + gtCount;
    }

    return count;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls
   * the buildKdTreePresorted function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode*> wherein each KdNode contains a (x,y,z,w...) tuple
   * numThreads - the number of threads
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
public:
  static KdNode<K,V,N>* createKdTree(std::vector<KdNode<K,V,N>*>& kdNodes,
                                     signed_size_t numThreads, signed_size_t maximumSubmitDepth) {

    struct timespec startTime, endTime;

    // Create a temporary vector for use in sorting the kdNodes vector and building the k-d tree.
    startTime = getTime();
    std::vector<KdNode<K,V,N>*> temporary(kdNodes.size());

    // Sort the kdNodes vector using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
    startTime = getTime();
    mergeSortReferenceAscending(kdNodes.data(), temporary.data(), 0, kdNodes.size() - 1,
                                0, maximumSubmitDepth, 0);
    endTime = getTime();
    double sortTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    // Remove KdNodes with duplicate tuples via one pass through the kdNodes vector.
    startTime = getTime();
    signed_size_t end = removeDuplicates(kdNodes.data(), 0, kdNodes.size());
    endTime = getTime();
    double removeTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    // Start the timer to time building the k-d tree.
    startTime = getTime();

    // Determine the maximum depth of the k-d tree, which is log2( kdNodes.size() ).
    signed_size_t size = kdNodes.size();
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the buildKdTree function because that coordinate depends only on the depth of
    // recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // Because the partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // (for e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference array
    // above using 0 as the leading key of the super key.
    std::vector<signed_size_t> permutation(maxDepth);
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      permutation[i] = i % K_DIMENSIONALITY;
    }

    // Build the k-d tree with multiple threads if possible.
    KdNode<K,V,N>* root = buildKdTreePresorted(kdNodes.data(), temporary.data(), permutation, 0, end,
                                               kdNodes.size(), maximumSubmitDepth);
    endTime = getTime();
    double kdTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    // Verify the k-d tree and report the number of kdNodes.
    startTime = getTime();
    signed_size_t numberOfNodes;
    numberOfNodes = root->verifyKdTree(maximumSubmitDepth, 0);
    endTime = getTime();
    double verifyTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
    std::cout << "Number of nodes = " << numberOfNodes << std::endl << std::endl;

    std::cout << "totalTime = " << std::fixed << std::setprecision(2) << (sortTime + removeTime + kdTime + verifyTime)
              << "  sortTime = " << sortTime << "  removeTime = " << removeTime
              << "  kdTime = " << kdTime << "  verifyTime = " << verifyTime << std::endl << std::endl;

    // Return the pointer to the root of the k-d tree.
    return root;
  }

  /*
   * Walk the k-d tree to delete each KdNode.
   */
public:
  void deleteKdTree() {

    // Delete the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->deleteKdTree();
    }
    // Delete the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->deleteKdTree();
    }
    // Delete the current KdNode.
    delete this;
  }

  /*
   * The insideBounds function determines whether KdNode::tuple lies inside the
   * hyper-rectangle defined by the query lower and upper bound vectors.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *
   * return true if inside, false if outside
   */
private:
  bool insideBounds(std::vector<K> const& queryLower, std::vector<K> const& queryUpper,
                    std::vector<bool> const& enable) {
    bool inside = true;
    for (signed_size_t i = 0; i < queryLower.size(); ++i) {
      if (enable[i] && (queryLower[i] > tuple[i] || queryUpper[i] < tuple[i])) {
        inside = false;
        break;
      }
    }
    return inside;
  }

  /*
   * The regionSearch function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * permutation - vector that specifies permutation of the partition coordinate
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   * enable - a vector that specifies the dimensions on which to prune the region search
   *
   * return a list that contains the KdNodes that lie within the cutoff distance of the query node
   */
private:
  std::list<KdNode<K,V,N>*> regionSearch(std::vector<K> const& queryLower, std::vector<K> const& queryUpper,
                                         std::vector<signed_size_t> const& permutation,
                                         signed_size_t maximumSubmitDepth, signed_size_t depth,
                                         std::vector<bool> const& enable) {

    // The partition cycles as x, y, z, w...
    signed_size_t p = permutation[depth];

    // If the KdNode is within the query hyper-rectangle for each of the k dimensions,
    // add the KdNode to the list of KdNodes that lie inside the hyper-cube. The
    // following loop is equivalent to the IN_REGION pseudo-Algol code proposed
    // by Jon Bentley in his CACM article.
    std::list<KdNode<K,V,N>*> result;
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_back(this);
    }

    // Determine whether to search the < and > branches of the k-d tree. Although
    // the superKeyCompare function can produce a different result for the == case
    // than does comparing only the leading keys of the super-keys, that result
    // will avoid unnecessary searching of a sub-tree (at the expense of a more
    // precise super-key comparison) but the unnecessary search/ appears not to
    // change the outcome of this recursive regionSearch function.
#ifdef NO_SUPER_KEY
    bool searchLT = ltChild != nullptr && (tuple[p] >= queryLower[p] || !enable[p]);
    bool searchGT = gtChild != nullptr && (tuple[p] <= queryUpper[p] || !enable[p]);;
#else
    bool searchLT = ltChild != nullptr && (superKeyCompare(tuple, queryLower.data(), p) >= 0 || !enable[p]);
    bool searchGT = gtChild != nullptr && (superKeyCompare(tuple, queryUpper.data(), p) <= 0 || !enable[p]);
#endif

    // Do both branches require searching and is a child thread available?
    if (searchLT && searchGT && maximumSubmitDepth >= 0 && depth <= maximumSubmitDepth) {

      // Yes, both branches of the tree require searching and a child thread is available,
      // so prepare to search the < branch with a child thread.
      std::future< std::list<KdNode<K,V,N>*> > searchFuture;

      // Search the < branch?
      if (searchLT) {
        
        // Yes, search the < branch asynchronously with a child thread.
        searchFuture =
          std::async(std::launch::async, [&] {
                                           return ltChild->regionSearch(queryLower, queryUpper, permutation,
                                                                        maximumSubmitDepth, depth + 1, enable);
                                         });

        // Search the > branch?
        std::list<KdNode<K,V,N>*> gtResult;
        if (searchGT) {
          
          // Yes, search the > branch  with the master thread.
          gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Get the result of searching the < branch with the child thread.
        std::list<KdNode<K,V,N>*> ltResult;
        try {
          ltResult = searchFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Append the results of searching the < and > branches to the result (if any) for this KdNode.
        result.splice(result.end(), ltResult);
        result.splice(result.end(), gtResult);

      } else {

        // No, don't search the < branch. Search the > branch?
        std::list<KdNode<K,V,N>*> gtResult;
        if (searchGT) {
          
          // Yes, search the > branch  with the master thread.
          gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Append the result of searching the > branch to the result (if any) for this KdNode.
        result.splice(result.end(), gtResult);
      }

    } else {
      
      // No, both branches do not require searching. Search the < branch with the master thread?
      if (searchLT) {
        auto ltResult = ltChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), ltResult); // Can't substitute regionSearch(...) for ltResult.
      }

      // Search the > branch with the master thread?
      if (searchGT) {
        auto gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), gtResult); // Can't substitute regionSearch(...) for gtResult.
      }

    }

    return result;
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * size - the number of points in the coordinates vector after removal of duplicates
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  std::list<KdNode<K,V,N>*> searchRegion(std::vector<K>& queryLower, std::vector<K>& queryUpper,
                                         signed_size_t maximumSubmitDepth, signed_size_t size) {
    
    // Determine the maximum depth of the k-d tree, which is log2(size).
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    std::vector<signed_size_t> permutation(maxDepth);
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % queryLower.size();
    }

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (signed_size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over all dimensions and return the resulting list of KdNodes.
    std::vector<bool> enable(queryLower.size(), true);
    return regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector

   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * size - the number of points in the coordinates vector after removal of duplicates
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  std::list<KdNode<K,V,N>*> searchRegion(std::vector<K>& queryLower, std::vector<K>& queryUpper,
                                         signed_size_t maximumSubmitDepth, signed_size_t size,
                                         std::vector<bool> const& enable) {
    
    // Determine the maximum depth of the k-d tree, which is log2(size).
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    std::vector<signed_size_t> permutation(maxDepth);
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % queryLower.size();
    }

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (signed_size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over the enabled dimensions and return the resulting list of KdNodes.
    return regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
  }

  /*
   * Walk the k-d tree and append to a list each KdNode that lies inside
   * the hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   *
   * return a list of KdNodes that lie within the query hyper-rectangle.
   */
public:
  std::list<KdNode<K,V,N>*> bruteRegion(std::vector<K> const& queryLower, std::vector<K> const& queryUpper) {

    // Append the KdNode to the list if it lies inside the query bounds.
    std::list<KdNode<K,V,N>*> result;
    std::vector<bool> enable(queryLower.size(), true);
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_back(this);
    }
    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      result.splice(result.end(), ltChild->bruteRegion(queryLower, queryUpper));
    }
    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      result.splice(result.end(), gtChild->bruteRegion(queryLower, queryUpper));
    }

    return result;
  }

  /*
   * Search the k-d tree for all possible M nearest geometric neighbors by adding them
   * to the NearestNeighborHeap.  Exclude a branch of the tree wherein it is guaranteed
   * that all nodes in that branch are farther way than the current farthest node stored
   * in the heap.
   *
   * Calling parameters:
   *
   * heap - an instance of NearestNeighborHeap
   * permutation - vector that specifies permutation of the partition coordinate
   * depth - depth in the k-d tree
   */
private:
  void nearestNeighbors(NearestNeighborHeap<K,V,N>& heap, std::vector<signed_size_t> const& permutation,
                        signed_size_t depth) {

    // The partition permutes as x, y, z, w...
    signed_size_t p = permutation.at(depth);

    // If query[p] < tuple[p], descend the < branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    if (heap.query[p] < tuple[p]) {
      if (ltChild != nullptr) {  // If not at the bottom of the tree, descend the < branch unconditionally.
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the > branch and then attempt to add the node to the heap.
      double dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (gtChild != nullptr) { // If not at the bottom of the tree, descend the > branch
          gtChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current KdNode to the heap.
      }
    }
    // If query[p] > tuple[p], descend the > branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    else if (heap.query[p] > tuple[p]) {
      if (gtChild != nullptr) {  // If not at the bottom of the tree, descend the > branch unconditionally.
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the < branch and then attempt to add the node to the heap.
      double dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (ltChild != nullptr) {
          ltChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current node to the heap.
      }
    }
    // Because query[p] == tuple[p], the probability of finding nearest neighbors is equal for both branches
    // of the tree, so descend both branches and then attempt to add the current node to the heap.
    else {
      if (ltChild != nullptr) {
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      if (gtChild != nullptr) {
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      heap.add(this);  // Attempt to add the current node to the heap.
    }
  }

  /*
   * Find M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   * size - the number of points in the coordinates vector after removal of duplicates
   *
   */
public:
  std::list< std::pair<double, KdNode<K,V,N>*> > findNearestNeighbors(std::vector<K> const& query,
                                                                      signed_size_t numNeighbors,
                                                                      signed_size_t size) {
    
    // Determine the maximum depth of the k-d tree, which is log2(size).
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    signed_size_t numDimensions = query.size();
    std::vector<signed_size_t> permutation(maxDepth);
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % numDimensions;
    }

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K,V,N> heap(query, numNeighbors);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K,V,N>*> > result;
    for (signed_size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }
  
  /*
   * Find M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   * size - the number of points in the coordinates vector after removal of duplicates
   * enable - a vector that specifies the dimensions on which to test distance
   *
   */
public:
  std::list< std::pair<double, KdNode<K,V,N>*> > findNearestNeighbors(std::vector<K> const& query, signed_size_t numNeighbors,
                                                                      signed_size_t size, std::vector<bool> const& enable) {
    
    // Determine the maximum depth of the k-d tree, which is log2(size).
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    signed_size_t numDimensions = query.size();
    std::vector<signed_size_t> permutation(maxDepth);
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % numDimensions;
    }

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K,V,N> heap(query, numNeighbors, enable);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K,V,N>*> > result;
    for (signed_size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }
  
  /*
   * Walk the k-d tree and attempt to add each KdNode to the NearestNeighborHeap.
   *
   * Calling parameter:
   *
   * heap - an instance of NearestNeighborHeap
   */
private:
  void allNeighbors(NearestNeighborHeap<K,V,N>& heap) {

    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->allNeighbors(heap);
    }
    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->allNeighbors(heap);
    }
    // Attempt to add the current KdNode to the heap.
    heap.add(this);
  }

  /*
   * Find M nearest neighbors to the query vector via brute force and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   *
   */
public:
  std::list< std::pair<double, KdNode<K,V,N>*> > bruteNearestNeighbors(std::vector<K> const& query,
                                                                       signed_size_t numNeighbors) {
    
    // Create the heap, walk the k-d tree, and attempt to add each KdNode to the heap.
    NearestNeighborHeap<K,V,N> heap(query, numNeighbors);
    allNeighbors(heap);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K,V,N>*> > result;
    for (signed_size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameters:
   *
   * tuple - the tuple to print
   * dim - the number of dimensions
   */
public:
  static void printTuple(K const* tuple)
    {
      std::cout << "(" << tuple[0] << ",";
      for (signed_size_t i = 1; i < N - 1; ++i) {
        std::cout << tuple[i] << ",";
        std::cout << tuple[N - 1] << ")";
      }
    }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameter:
   *
   * tuple - the tuple as a vector
   */
public:
  static void printTuple(std::vector<K> const& tuple)
    {
      std::cout << "(" << tuple[0] << ",";
      for (signed_size_t i = 1; i < tuple.size() - 1; ++i) {
        std::cout << tuple[i] << ",";
        std::cout << tuple[tuple.size() - 1] << ")";
      }
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
  void printKdTree(signed_size_t depth)
    {
      if (gtChild != nullptr) {
        gtChild->printKdTree(depth + 1);
      }
      for (signed_size_t i = 0; i < depth; ++i) std::cout << "       ";
      printTuple(tuple);
      std::cout << std::endl;
      if (ltChild != nullptr) {
        ltChild->printKdTree(depth + 1);
      }
    }
}; // class KdNode

/*
 * The NearestNeighborHeap Class implements a fixed length heap of both containing both a KdNode and euclidean distance
 * from the tuple in the node to a query point.  When a KdNode is added to the heap it is unconditionally placed in
 * the heap until the heap is full.  After the heap is full, a KdNode is added to the heap only if the calculated
 * distance from the query point to the tuple is less than the farthest KdNode currently in the heap; and in that
 * case, the current farthest KdNode and distance are removed from the heap to make room for it.
 *
 * The heap is maintained in two corresponding fixed length vectors, one for the KdNodes and one for the distance to
 * the query point.  These vectors are stored in order of decreasing distance.  When a new node is added, regardless
 * of whether or not the heap is full, a heap sort is done to place the new KdNode in the proper order in the heap.
 * Hence, the farthest KdNode is always at the top of the heap (index 1).
 *
 * For a discussion of heap sort and a priority queue implemented via a heap, see Section 2.4 "Priority Queues"
 * pp. 308-335 in "Algorithms Fourth Edition" by Robert Sedgewick and Kevin Wayne, Addison-Wesley, 2011.
 */
template <typename K, typename V, signed_size_t N>
class NearestNeighborHeap {
public:
  std::vector<K> query; // query point for which nearest neighbors will be found
  std::vector<bool> enable;
private:
  signed_size_t reqDepth; // requested number of nearest neighbors
  std::vector<KdNode<K,V,N>* > nodes; // vector of KdNodes that are the nearest neighbors
  std::vector<double> dists; // vector of squared distances
  signed_size_t curDepth; // number of nearest nodes/distances on the heap

  /*
   * Constructor that enables distance test on all dimensions
   *
   * Calling parameters:
   *
   * query - a vector that defines the query point
   * numNeighbors - the number of nearest neighbors desired
   */
public:
  NearestNeighborHeap(std::vector<K> const& query, signed_size_t numNeighbors) {
    this->nodes.resize(numNeighbors + 1, nullptr); // heap of KdNode* (address 0 is unused)
    this->dists.resize(numNeighbors + 1, 0); // corresponding heap of distances (initialized to 0)
    this->reqDepth = numNeighbors;
    this->curDepth = 0;
    this->query = query;
    this->enable.assign(query.size(), true);
  }
  
  /*
   * Constructor that enables distance test on only specified dimensions
   *
   * Calling parameters:
   *
   * query - a vector that defines the query point
   * numNeighbors - the number of nearest neighbors desired
   * enable - a vector that specifies the dimensions on which to test distance
   */
public:
  NearestNeighborHeap(std::vector<K> const& query, signed_size_t numNeighbors, std::vector<bool> const& enable) {
    this->nodes.resize(numNeighbors + 1, nullptr); // heap of KdNode* (address 0 is unused)
    this->dists.resize(numNeighbors + 1, 0); // corresponding heap of distances (initialized to 0)
    this->reqDepth = numNeighbors;
    this->curDepth = 0;
    this->query = query;
    this->enable = enable;
  }

  /*
   * Swap two elements in the heap.
   *
   * Calling parameters:
   *
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  void swap(signed_size_t i, signed_size_t j) {
    double tempDist = dists[i];
    KdNode<K,V,N>* tempNode = nodes[i];
    dists[i] = dists[j];
    nodes[i] = nodes[j];
    dists[j] = tempDist;
    nodes[j] = tempNode;
  }

  /*
   * Allow an element to rise upward through the heap.
   *
   * Calling parameter:
   *
   * k - the index of the element
   */
private:
  void rise(signed_size_t k) {
    while (k > 1 && dists[k/2] < dists[k]) {
      swap(k/2, k);
      k = k/2;
    }
  }
  /*
   * Allow an element to fall downward through the heap.
   *
   * Calling parameter:
   *
   * k - the index of the element
   */
private:
  void fall(signed_size_t k) {
    while (2*k <= curDepth) {
      int j = 2*k;
      if (j < curDepth && dists[j] < dists[j+1]) {
        ++j;
      }
      if (dists[k] >= dists[j]) {
        break;
      }
      swap(k, j);
      k = j;
    }
  }

  /*
   * Remove the top element of the heap and re-order the remaining elements.
   *
   * return a pair that contains a pointer to the top KdNode and the distance to that KdNode
   */
public:
  std::pair<double, KdNode<K,V,N>*> removeTop() {
    std::pair<double, KdNode<K,V,N>*> returnPair = std::make_pair(dists[1], nodes[1]);
    swap(1, curDepth--);
    nodes[curDepth+1] = nullptr;
    fall(1);
    return returnPair;
  }
  
  /*
   * Add a new KdNode to the NearestNeighborHeap if its distance to the
   * query point is less than curMaxDistance.
   *
   * Calling parameter:
   *
   * node - KdNode to potentially add to the heap
   */
public:
  void add(KdNode<K,V,N>* node) {
    // Find the distance by subtracting the query from the tuple and
    // calculating the sum of the squared distances. Note that conversion
    // from type K to double may result in loss of precision but avoids
    // the possibility of integer overflow.
    double dist2 = 0.0;
    for (int i = 0; i < query.size(); ++i) {
      // Add the squared coordinate distance only if the dimension is enabled.
      if (enable[i]) {
        K comp = node->tuple[i] - query[i];
        dist2 += static_cast<double>(comp) * static_cast<double>(comp);
      }
    }
    // If the queue is not full, add the point to the bottom of the heap unconditionally and let it rise;
    if (!heapFull()) {
      dists[++curDepth] = dist2;
      nodes[curDepth] = node;
      rise(curDepth);
    }
    // otherwise, if the point is closer than the top of the heap, overwrite the top and let it fall.
    else if (dist2 < curMaxDist()) {
      dists[1] = dist2;
      nodes[1] = node;
      fall(1);
    }
    return;
  }

  /* Return the current maximum distance, i.e., dists[1] */
public:
  double curMaxDist() {
    return dists[1];
  }

  /* Return true if the heap is full. */
public:
  bool heapFull() {
    return curDepth >= reqDepth;
  }
}; // class NearestNeighborHeap

/*
 * The randomLongInInterval function creates a random long in the interval [min, max].  See
 * http://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
 *
 * Calling parameters:
 *
 * min - the minimum long value desired
 * max - the maximum long value desired
 *
 * returns: a random long
 */
static kdKey_t randomLongInInterval(kdKey_t min, kdKey_t max) {
  // subtract 32768 from range to avoid overflows.
  return min + (kdKey_t)((((double)rand()) / ((double)RAND_MAX)) * (max - min - 32768));
}


#ifdef TEST_KD_TREE
/* Create a simple k-d tree and print its topology for inspection. */
int main(int argc, char** argv)
{
  struct timespec startTime, endTime;

  // Set the defaults then parse the input arguments.
  signed_size_t numPoints = 262144;
  signed_size_t numNeighbors = 5;
  signed_size_t extraPoints = 100;
  signed_size_t numThreads = 4;
  bool bruteForceSearch = false;
  bool bruteForceRegion = false;
  signed_size_t maximumNumberOfNodesToPrint = 5;
  kdKey_t searchDistance = 1000000000000000000L;

  for (signed_size_t i = 1; i < argc; ++i) {
    if (0 == strcmp(argv[i], "-n") || 0 == strcmp(argv[i], "--numPoints")) {
      numPoints = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-m") || 0 == strcmp(argv[i], "--numNeighbors")) {
      numNeighbors = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-x") || 0 == strcmp(argv[i], "--extraPoints")) {
      extraPoints = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-t") || 0 == strcmp(argv[i], "--numThreads")) {
      numThreads = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-s") || 0 == strcmp(argv[i], "--searchDistance")) {

      continue;
    }
    if (0 == strcmp(argv[i], "-p") || 0 == strcmp(argv[i], "--maximumNodesToPrint")) {
      maximumNumberOfNodesToPrint = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-b") || 0 == strcmp(argv[i], "--bruteForceSearch")) {
      bruteForceSearch = !bruteForceSearch;
      continue;
    }
    if (0 == strcmp(argv[i], "-c") || 0 == strcmp(argv[i], "--bruteForceRegion")) {
      bruteForceRegion = !bruteForceRegion;
      continue;
    }
    std::cout << "illegal command-line argument: " << argv[i] << std::endl;
    exit(1);
  }

  // Declare and initialize the kdNodes vector and initialize it with tuples,
  // for example (x,y,z,w), in the half-open interval [0, MAX] where MAX is the
  // maximum value for the kdKey_t type. Create extraPoints-1 duplicate elements,
  // where extraPoints <= numPoints, to test the removal of duplicate points.
  extraPoints = (extraPoints <= numPoints) ? extraPoints : numPoints;
  std::vector<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> kdNodes(numPoints + extraPoints - 1);
  // Create each KdNode with a tuple (all elements == 0) and a value.
  for (signed_size_t i = 0; i < kdNodes.size(); ++i) {
    kdNodes[i] = new KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>(i); // KdNode::value is its vector index.
  }
  // Initialize each KdNode::tuple except for the extra points.
  for (signed_size_t i = 0; i < numPoints; ++i) {
    for (signed_size_t j = 0; j < K_DIMENSIONALITY; ++j) {
      kdNodes[i]->tuple[j] = randomLongInInterval(0, std::numeric_limits<kdKey_t>::max());
    }
  }
  // Replicate the tuples from the uppermost KdNodes to obtain the extra KdNodes' duplicate tuples.
  for (signed_size_t i = 1; i < extraPoints; ++i) {
    for (signed_size_t j = 0; j < K_DIMENSIONALITY; ++j) {
      kdNodes[numPoints - 1 + i]->tuple[j] = kdNodes[numPoints - 1 - i]->tuple[j];
    }
  }

  // Calculate the number of child threads to be the number of threads minus 1, then
  // calculate the maximum tree depth at which to launch a child thread.  Truncate
  // this depth such that the total number of threads, including the master thread, is
  // an integer power of 2, hence simplifying the launching of child threads by restricting
  // them to only the < branch of the tree for some depth in the tree.
  signed_size_t n = 0;
  if (numThreads > 0) {
    while (numThreads > 0) {
      ++n;
      numThreads >>= 1;
    }
    numThreads = 1 << (n - 1);
  }
  else {
    numThreads = 0;
  }
  signed_size_t childThreads = numThreads - 1;
  signed_size_t maximumSubmitDepth = -1;
  if (numThreads < 2) {
    maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
  }
  else if (numThreads == 2) {
    maximumSubmitDepth = 0;
  }
  else {
    maximumSubmitDepth = static_cast<signed_size_t>(floor(log(static_cast<double>(childThreads)) / log(2.)));
  }
  std::cout << std::endl << "Max number of threads = " << numThreads << "  max submit depth = "
            << maximumSubmitDepth << std::endl << std::endl;

  // Create the k-d tree.
  auto root = KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>::createKdTree(kdNodes, numThreads, maximumSubmitDepth);

  // Search the k-d tree via region search for the KdNodes that lie within a hyper-cube centered near the origin.
  std::vector<kdKey_t> query(K_DIMENSIONALITY);
  std::vector<kdKey_t> queryLower(K_DIMENSIONALITY);
  std::vector<kdKey_t> queryUpper(K_DIMENSIONALITY);
  for (signed_size_t i = 0; i < K_DIMENSIONALITY; ++i) {
    query[i] = i;
    queryLower[i] = query[i] - searchDistance;
    queryUpper[i] = query[i] + searchDistance;
  }
  startTime = getTime();
  auto regionFast = root->searchRegion(queryLower, queryUpper, maximumSubmitDepth, kdNodes.size());
  endTime = getTime();
  double fastRegionTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  std::cout << "fast region time = " << std::fixed << std::setprecision(6) << fastRegionTime << " seconds" << std::endl << std::endl;

  std::cout << regionFast.size() << " nodes within " << searchDistance << " units of ";
  KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>::printTuple(query);
  std::cout << " in all dimensions." << std::endl << std::endl;
  if (regionFast.size() != 0) {
    regionFast.sort();
    signed_size_t maxNodesToPrint = maximumNumberOfNodesToPrint;
    std::cout << "List of the first <= " << maximumNumberOfNodesToPrint << " fast k-d nodes within a "
              << searchDistance << "-unit search distance follows:" << std::endl << std::endl;
    for (std::list<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*>::iterator it = regionFast.begin(); it != regionFast.end(); ++it) {
      KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>::printTuple((*it)->getTuple());
      std::cout << std::endl;
      --maxNodesToPrint;
      if (maxNodesToPrint == 0) {
        break;
      }
    }
    std::cout << std::endl;
  }

  // Verify that no duplicate KdNodes exist on the list returned from region search.
  auto itr1 = regionFast.begin();
  auto itr2 = itr1;
  ++itr2;
  for ( ; itr2 != regionFast.end(); ++itr1, ++itr2) {
    if (*itr1 == *itr2) {
      throw std::runtime_error("Duplicate KdNode* on region-search list\n");
    }
  }

  // Search the k-d tree via brute force for the KdNodes that lie within a hyper-cube centered near the origin.
  if (bruteForceRegion) {
    startTime = getTime();
    auto regionSlow = root->bruteRegion(queryLower, queryUpper);
    endTime = getTime();
    double slowRegionTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    std::cout << "slow region time = " << std::fixed << std::setprecision(6) << slowRegionTime << " seconds" << std::endl << std::endl;

    std::cout << regionSlow.size() << " nodes within " << searchDistance << " units of ";
    KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>::printTuple(query);
    std::cout << " in all dimensions." << std::endl << std::endl;
    if (regionSlow.size() != 0) {
      regionSlow.sort();
      signed_size_t maxNodesToPrint = maximumNumberOfNodesToPrint;
      std::cout << "List of the first <= " << maximumNumberOfNodesToPrint << " slow k-d nodes within a "
                << searchDistance << "-unit search distance follows:" << std::endl << std::endl;
      for (std::list<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*>::iterator it = regionSlow.begin(); it != regionSlow.end(); ++it) {
        KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>::printTuple((*it)->getTuple());
        std::cout << std::endl;
        --maxNodesToPrint;
        if (maxNodesToPrint == 0) {
          break;
        }
      }
      std::cout << std::endl;
    }

    // Verify that the region-search and brute-force lists are identical.
    auto itrf = regionFast.begin();
    for (auto itrs = regionSlow.begin(); itrs != regionSlow.end(); ++itrf, ++itrs) {
      if (*itrf != *itrs) {
        throw std::runtime_error("Non-identical region-search and brute-force lists\n");
      }
    }
  }

  // It is impossible to find more nearest neighbors than there are points.
  numNeighbors = std::min(numNeighbors, numPoints + extraPoints + 1);

  // Search the k-d tree for the numNeighbors nearest neighbors to the first tuple.
  startTime = getTime();
  auto neighborsFast = root->findNearestNeighbors(query, numNeighbors, kdNodes.size());
  endTime = getTime();
  double fastNeighborTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  std::cout << "fast neighbor time = " << std::fixed << std::setprecision(6) << fastNeighborTime << " seconds" << std::endl << std::endl;
  std::cout << "fast neighbor list size = " << neighborsFast.size() << std::endl << std::endl;

  // Find nearest neighbors via brute force if requested.
  if (bruteForceSearch) {
    startTime = getTime();
    auto neighborsSlow = root->bruteNearestNeighbors(query, numNeighbors);
    endTime = getTime();
    double slowNeighborTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    std::cout << "slow neighbor time = " << std::fixed << std::setprecision(6) << slowNeighborTime << " seconds" << std::endl << std::endl;
    std::cout << "slow neighbor list size = " << neighborsSlow.size() << std::endl << std::endl;

    auto itf1 = neighborsFast.begin();
    auto its1 = neighborsSlow.begin();
    
    // Compare the first k-d tree distance to the first brute-force distance.
    if (itf1->first != its1->first) {
      char msg[256];
      snprintf(msg, 256, "fast distance[0] = %f  !=  slow distance[0] = %f\n", itf1->first, its1->first);
      throw std::runtime_error(msg);
    }
    
    // Compare the first k-d tree KdNode* to the first brute-force KdNode*.
    if (itf1->second != its1->second) {
      char msg[256];
      snprintf(msg, 256, "fast KdNode*[0] != slow KdNode*[0]\n");
      throw std::runtime_error(msg);
    }

    std::cout << "fast and slow distances squared follow:" << std::endl << std::endl;
    std::cout << std::fixed << std::setprecision(0) << itf1->first << "\t" << its1->first << std::endl;

    // Compare the remaining distances and KdNode*
    auto itf2 = itf1;
    auto its2 = its1;
    ++itf2;
    ++its2;
    signed_size_t i = 1;
    for ( ; itf2 != neighborsFast.end(); ++itf1, ++its1, ++itf2, ++its2, ++i) {
      // Print the fast and slow distances for this iteration of the loop.
      std::cout << std::fixed << std::setprecision(0) << itf2->first << "\t" << its2->first << std::endl;
      // Ensure that the fast distances increase monotonically.
      if (itf1->first > itf2->first) {
        char msg[256];
        snprintf(msg, 256, "fast distance[%ld] = %f  >  fast distance[%ld] = %f\n",
                 i - 1, itf1->first, i, itf2->first);
        throw std::runtime_error(msg);
      }
      // Ensure that the slow distances increase monotonically.
      if (its1->first > its2->first) {
        char msg[256];
        snprintf(msg, 256, "slow distance[%ld] = %f  >  slow distance[%ld] = %f\n",
                 i-1, its1->first, i, its2->first);
        throw std::runtime_error(msg);
      }
      // Compare the ith k-d tree distance to the ith brute-force distance.
      if (itf2->first != its2->first) {
        char msg[256];
        snprintf(msg, 256, "fast distance[%ld] = %f  !=  slow distance[%ld] = %f\n",
                 i, itf2->first, i, its2->first);
        throw std::runtime_error(msg);
      }
      // Compare the ith k-d tree KdNode* to the ith brute-force KdNode*.
      if (itf2->second != its2->second) {
        char msg[256];
        snprintf(msg, 256, "fast KdNode*[%ld] != slow KdNode*[%ld]\n", i, i);
        throw std::runtime_error(msg);
      }
    }
    std::cout << std::endl;
  }

  // Delete the k-d tree.
  startTime = getTime();
  root->deleteKdTree();
  endTime = getTime();
  double deleteTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  std::cout << "deleteTime = " << std::fixed << std::setprecision(6) << deleteTime << " seconds" << std::endl << std::endl;

  return 0;
}
#endif // #ifdef TEST_KD_TREE
