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
 * n elements of data, a balanced k-d tree is built in O(kn log n) + O((k-1)n log n)
 * time by first sorting the data in each of k dimensions, then building the k-d tree
 * in a manner that preserves the order of the k sorts while recursively partitioning
 * the data at each level of the k-d tree. No further sorting is necessary.
 *
 * Gnu g++ compilation options are: -lm -O3 -std=c++11 -DTEST_KD_TREE
 *
 * Optional compilation options are:
 *
 * -D K_DIMENSIONALITY=n - The k-dimensionality of the k-d tree (default 3).
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 * -D MACH - Use a Mach equivalent to the clock_gettime(CLOCK_REALTIME, &time) function
 *           but this option appears to no longer be necessary.
 */

#include <limits>
#include <math.h>
#include <time.h>
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
      K compare = superKeyCompare(kdNodes[j]->tuple, kdNodes[j - 1]->tuple, i);
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
  static KdNode<K,V,N>* buildKdTree(KdNode<K,V,N>*** references,
                                    std::vector< std::vector<signed_size_t> > const& permutation,
                                    signed_size_t start, signed_size_t end,
                                    signed_size_t maximumSubmitDepth, signed_size_t depth) {

    KdNode<K,V,N>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    signed_size_t p = permutation.at(depth).at(permutation.at(0).size() - 1);

    // Obtain the reference array that corresponds to the most significant key.
    KdNode<K,V,N>** reference = references[permutation.at(depth).at(N)];

    if (end == start) {

      // Only one reference was passed to this function, so add it to the tree.
      node = reference[end];

    }
    else if (end == start + 1) {

      // Two references were passed to this function in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child.
      node = reference[start];
      node->gtChild = reference[end];

    }
    else if (end == start + 2) {

      // Three references were passed to this function in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = reference[start + 1];
      node->ltChild = reference[start];
      node->gtChild = reference[end];

    }
    else if (end > start + 2) {

      // Four or more references were passed to this function, so the
      // median element of the reference array is chosen as the tuple
      // about which the other reference arrays will be partitioned
      // Avoid overflow when computing the median.
      signed_size_t median = start + ((end - start) / 2);

      // Store the median element of the reference array in a new KdNode.
      node = reference[median];

      // Build both branches with child threads at as many levels of the tree
      // as possible.  Create the child threads as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so one thread will be used.
        // Initialize startIndex=1 so that the 'for' loop that partitions the
        // reference arrays will partition a number of arrays equal to N.
        signed_size_t startIndex = 1;

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNode<K,V,N>** dst = references[permutation.at(depth).at(0)];
          KdNode<K,V,N>** tmp = references[permutation.at(depth).at(1)];
          for (int i = start; i <= end; ++i) {
            dst[i] = reference[i];
          }
          // Sort the lower half of references[permut[0]] with the current thread.
          mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1, maximumSubmitDepth, depth);
          // Sort the upper half of references[permut[0]] with the current thread.
          mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1, maximumSubmitDepth, depth);
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        K* tuple = node->tuple;
        for (signed_size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNode<K,V,N>** src = references[permutation.at(depth).at(i)];
          KdNode<K,V,N>** dst = references[permutation.at(depth).at(i - 1)];

          // Fill the lower and upper halves of one reference array
          // in ascending order with the current thread.
          for (signed_size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
            KdNode<K,V,N>* src_j = src[j];
            K compare = superKeyCompare(src_j->tuple, tuple, p);
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
        // reference arrays will partition a number of arrays equal to N.
        signed_size_t startIndex = 1;

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNode<K,V,N>** dst = references[permutation.at(depth).at(0)];
          KdNode<K,V,N>** tmp = references[permutation.at(depth).at(1)];
          // Copy and sort the lower half of references[permut[0]] with a child thread.
          std::future<void> copyFuture =
            std::async(std::launch::async, [&] {
                                   for (int i = start; i <= median - 1; ++i) {
                                     dst[i] = reference[i];
                                   }
                                   mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1, maximumSubmitDepth, depth);
                                 });

          // Copy and sort the upper half of references[permut[0]] with the current thread.
          for (int i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
          }
          mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1, maximumSubmitDepth, depth);

          // Wait for the child thread to finish execution.
          try {
            copyFuture.get();
          }
          catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
        }

        // Create a copy of the node->tuple array so that the current thread
        // and the child thread do not contend for read access to this array.
        K* tuple = node->tuple;
        K* point = new K[N];
        for (signed_size_t i = 0; i < N; ++i) {
          point[i] = tuple[i];
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        for (signed_size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNode<K,V,N>** src = references[permutation.at(depth).at(i)];
          KdNode<K,V,N>** dst = references[permutation.at(depth).at(i - 1)];

          // Two threads may be used to partition the reference arrays, analogous to
          // the use of two threads to merge the results for the merge sort algorithm.
          // Fill one reference array in ascending order with a child thread.
          std::future<void> partitionFuture =
            std::async(std::launch::async, [&] {
                                   for (signed_size_t lower = start - 1, upper = median, j = start; j <= median; ++j) {
                                     KdNode<K,V,N>* src_j = src[j];
                                     K compare = superKeyCompare(src_j->tuple, point, p);
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
            KdNode<K,V,N>* src_k = src[k];
            K compare = superKeyCompare(src_k->tuple, tuple, p);
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
          catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
        }

        // Delete the point array.
        delete[] point;

        // Recursively build the < branch of the tree with a child thread.
        // The recursive call to buildKdTree must be placed in a lambda
        // expression because buildKdTree is a template not a function.
        std::future<KdNode<K,V,N>*> buildFuture =
          std::async(std::launch::async, [&] {
                                 return buildKdTree(references, permutation, start, median - 1,
                                                    maximumSubmitDepth, depth + 1);
                               });

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, median + 1, end,
                                    maximumSubmitDepth, depth + 1);

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
  static void swap(std::vector<signed_size_t>& a, signed_size_t i, signed_size_t j) {
    signed_size_t t = a.at(i);
    a.at(i) = a.at(j);
    a.at(j) = t;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode*> wherein each KdNode contains a (x,y,z,w...) tuple
   * numDimensions - the number of dimensions
   * numThreads - the number of threads
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
public:
  static KdNode<K,V,N>* createKdTree(std::vector<KdNode<K,V,N>*>& kdNodes,
                                     signed_size_t numThreads, signed_size_t maximumSubmitDepth) {

    struct timespec startTime, endTime;

    // Create the references arrays including one additional array for use in building the k-d tree.
    KdNode<K,V,N>*** references = new KdNode<K,V,N>**[N + 1];

    // The first references array is the .data() array of the kdNodes vector.
    references[0] = kdNodes.data();

    // Allocate the remaining references arrays.
    for (int i = 1; i < N + 1; ++i) {
      references[i] = new KdNode<K,V,N>*[kdNodes.size()];
    }

    // Sort the first reference array using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
    startTime = getTime();
    mergeSortReferenceAscending(references[0], references[N], 0, kdNodes.size() - 1,
                                0, maximumSubmitDepth, 0);
    endTime = getTime();
    double sortTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    // Remove references to duplicate coordinates via one pass through the first reference array.
    startTime = getTime();
    signed_size_t end = removeDuplicates(references[0], 0, kdNodes.size());
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

    // It is unnecessary to compute either the permutation of the reference array or
    // the partition coordinate upon each recursive call of the buildKdTree function
    // because both depend only on the depth of recursion, so they may be pre-computed.
    // Create and initialize an 'indices' vector for the permutation calculation.
    // Because this vector is initialized with 0, 1, 2, 3, 0, 1, 2, 3, etc. (for
    // e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference
    // array above using 0 as the leading key of the super key.
    std::vector<signed_size_t> indices(N + 2);
    for (signed_size_t i = 0; i < indices.size() - 1; ++i) {
      indices[i] = i;
    }

    // Create a 'permutation' vector from the 'indices' vector to specify permutation
    // of the reference arrays and of the partition coordinate.
    std::vector< std::vector<signed_size_t> > permutation(maxDepth, std::vector<signed_size_t>(N + 2));

    // Fill the permutation vector by calculating the permutation of the indices vector
    // and the the partition coordinate of the tuple at each depth in the tree.
    for (signed_size_t i = 0; i < permutation.size(); ++i) {
      // The last entry of the indices vector contains the partition coordinate.
      indices.at(N + 1) = i % N;
      // Swap the first and second to the last elements of the indices vector.
      swap(indices, 0, N);
      // Copy the indices vector to one row of the permutation vector.
      permutation.at(i) = indices;
      // Swap the third and second to the last elements of the indices vector.
      swap(indices, N - 1, N);
    }

    // Build the k-d tree with multiple threads if possible.
    KdNode<K,V,N>* root = buildKdTree(references, permutation, 0, end,
                                      maximumSubmitDepth, 0);
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

    // Delete all but the first of the references arrays.
    for (int i = 1; i < N + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the root of the k-d tree.
    return root;
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
  std::list< std::pair<double, KdNode<K,V,N>*> > findNearestNeighbors(std::vector<K> const& query,
                                                            signed_size_t numNeighbors,
                                                            signed_size_t size,
                                                            std::vector<bool> const& enable) {
    
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
  std::list< std::pair<double, KdNode<K,V,N>*> > bruteNearestNeighbors(std::vector<K> const& query, signed_size_t numNeighbors) {
    
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
 * The NearestNeighborHeap Class implements a std::fixed length heap of both containing both a KdNode and euclidean distance
 * from the tuple in the node to a query point.  When a KdNode is added to the heap it is unconditionally placed in
 * the heap until the heap is full.  After the heap is full, a KdNode is added to the heap only if the calculated
 * distance from the query point to the tuple is less than the farthest KdNode currently in the heap; and in that
 * case, the current farthest KdNode and distance are removed from the heap to make room for it.
 *
 * The heap is maintained in two corresponding std::fixed length vectors, one for the KdNodes and one for the distance to
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
  }

  std::cout << std::endl;

  return 0;
}
#endif // #ifdef TEST_KD_TREE
