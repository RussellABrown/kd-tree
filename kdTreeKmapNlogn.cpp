/*
 * Copyright (c) 2015, 2021, 2023 Russell A. Brown
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
 * Gnu g++ compilation options are: -lm -O3 -pthread -D TEST_KD_TREE -W
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

#include <exception>
#include <forward_list>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>

using std::async;
using std::cout;
using std::endl;
using std::distance;
using std::exception;
using std::fixed;
using std::forward_list;
using std::future;
using std::launch;
using std::list;
using std::map;
using std::make_pair;
using std::min;
using std::numeric_limits;
using std::ostringstream;
using std::pair;
using std::ref;
using std::runtime_error;
using std::scientific;
using std::setprecision;
using std::streamsize;
using std::string;
using std::vector;

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
typedef streamsize signed_size_t;

/*
 * These are the types used for the test. Change the intrisic types in
 * these typedefs to test the k-d tree with different intrisic types.
 */
typedef int64_t kdKey_t;
typedef string kdValue_t;

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
  KdNode(V const& value) {
    
    this->value = value;
    this->ltChild = this->gtChild = this->duplicates = nullptr; // redundant
  }

public:
  ~KdNode() {
    // Delete each KdNode from the duplicates list.
    auto nextPtr = this->duplicates;
    while (nextPtr != nullptr) {
      auto tempPtr = nextPtr;
      nextPtr = nextPtr->duplicates;
      tempPtr->duplicates = nullptr; // Prevent recursive deletion.
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
   * dim - the number of dimensions
   *
   * returns a K result of comparing two K arrays
   */
private:
  inline
  static K superKeyCompare(K const* a,
                           K const* b,
                           signed_size_t const p,
                           signed_size_t const dim) {
    
    // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
    K diff = a[p] - b[p];
    for (signed_size_t i = 1; diff == 0 && i < dim; i++) {
      signed_size_t r = i + p;
      // A fast alternative to the modulus operator for (i + p) < 2 * dim.
      r = (r < dim) ? r : r - dim;
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
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceAscending(KdNode<K,V,N>** const reference,
                                          KdNode<K,V,N>** const temporary,
                                          signed_size_t const low,
                                          signed_size_t const high,
                                          signed_size_t const p,
                                          signed_size_t const dim,
                                          signed_size_t const maximumSubmitDepth,
                                          signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in ascending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) < 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in ascending order.
        auto sortFuture = async(launch::async, mergeSortTemporaryAscending, reference, temporary,
                                low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for sort future in mergeSortReferenceAscending\n");
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in ascending order.
        auto mergeFuture =
          async(launch::async, [&] {
                                 for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                   reference[k] =
                                     (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) <= 0)
                                     ? temporary[i++] : temporary[j--];
                                 }
                               });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in ascending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) > 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for merge future in mergeSortReferenceAscending\n");
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        auto const tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p, dim) > 0; --j) {
          reference[j] = reference[j - 1];
        }
        reference[j] = tmp;
      }
    }
  }

  /*
   * The mergeSortReferenceDescending function recursively subdivides the reference array then
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
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceDescending(KdNode<K,V,N>** const reference,
                                           KdNode<K,V,N>** const temporary,
                                           signed_size_t const low,
                                           signed_size_t const high,
                                           signed_size_t const p,
                                           signed_size_t const dim,
                                           signed_size_t const maximumSubmitDepth,
                                           signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in descending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) > 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in descending order.
        auto sortFuture = async(launch::async, mergeSortTemporaryDescending, reference, temporary,
                                low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for sort future in mergeSortReferenceDescending\n");
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in descending order.
        auto mergeFuture =
          async(launch::async, [&] {
                                 for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                   reference[k] =
                                     (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) >= 0)
                                     ? temporary[i++] : temporary[j--];
                                 }
                               });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in descending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p, dim) < 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for merge future in mergeSortReferenceDescending\n");
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        auto const tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p, dim) < 0; --j) {
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
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryAscending(KdNode<K,V,N>** const reference,
                                          KdNode<K,V,N>** const temporary,
                                          signed_size_t const low,
                                          signed_size_t const high,
                                          signed_size_t const p,
                                          signed_size_t const dim,
                                          signed_size_t const maximumSubmitDepth,
                                          signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in ascending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) < 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in ascending order.
        auto sortFuture = async(launch::async, mergeSortReferenceAscending, reference, temporary,
                                low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for sort future in mergeSortTemporaryAscending\n");
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in ascending order.
        auto mergeFuture =
          async(launch::async, [&] {
                                 for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                   temporary[k] =
                                     (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) <= 0)
                                     ? reference[i++] : reference[j--];
                                 }
                               });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in ascending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) > 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for merge future in mergeSortTemporaryAscending\n");
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
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p, dim) > 0) {
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
   * The mergeSortTemporaryDescending function recursively subdivides the reference array
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
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryDescending(KdNode<K,V,N>** const reference,
                                           KdNode<K,V,N>** const temporary,
                                           signed_size_t const low,
                                           signed_size_t const high,
                                           signed_size_t const p,
                                           signed_size_t const dim,
                                           signed_size_t const maximumSubmitDepth,
                                           signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in descending order.
        for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) > 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in descending order.
        auto sortFuture = async(launch::async, mergeSortReferenceDescending, reference, temporary,
                                low, mid, p, dim, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, dim, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for sort future in mergeSortTemporaryDescending\n");
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in descending order.
        auto mergeFuture =
          async(launch::async, [&] {
                                 for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                   temporary[k] =
                                     (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) >= 0)
                                     ? reference[i++] : reference[j--];
                                 }
                               });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in descending order.
        for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p, dim) < 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for merge future in mergeSortTemporaryDescending\n");
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
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p, dim) < 0) {
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
   * dim - the number of dimensions
   *
   * returns the end index of the reference array following removal of duplicate elements
   */
private:
  inline
  static signed_size_t removeDuplicates(KdNode<K,V,N>** const kdNodes,
                                        signed_size_t const i,
                                        signed_size_t const dim,
                                        signed_size_t const size) {
    
    signed_size_t end = 0;
    for (signed_size_t j = 1; j < size; ++j) {
      K const compare = superKeyCompare(kdNodes[j]->tuple, kdNodes[j - 1]->tuple, i, dim);
      if (compare < 0) {
        ostringstream buffer;
        buffer << "\n\nmerge sort failure: superKeyCompare(ref[" << j << "], ref["
               << end << "], " << i << ") = " << compare << "in removeDuplicates\n";
        throw runtime_error(buffer.str());
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
  static void swap(KdNode<K,V,N>** const a,
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
  static KdNode<K,V,N>* select_0_2(KdNode<K,V,N>* const a,
                                   KdNode<K,V,N>* const b,
                                   signed_size_t const p,
                                   signed_size_t const dim) {
    
    if (superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_2(KdNode<K,V,N>* const a,
                                   KdNode<K,V,N>* const b,
                                   signed_size_t const p,
                                   signed_size_t const dim) {
    
    if (superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_3_ab(KdNode<K,V,N>* const a,
                                      KdNode<K,V,N>* const b,
                                      KdNode<K,V,N>* const c,
                                      signed_size_t const p,
                                      signed_size_t const dim) {
    
    if (superKeyCompare(b->tuple, c->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_3(KdNode<K,V,N>* const a,
                                   KdNode<K,V,N>* const b,
                                   KdNode<K,V,N>* const c,
                                   signed_size_t const p,
                                   signed_size_t const dim) {
    
    if (superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_4_ab_cd(KdNode<K,V,N>* const a,
                                         KdNode<K,V,N>* const b,
                                         KdNode<K,V,N>* const c,
                                         KdNode<K,V,N>* const d,
                                         signed_size_t const p,
                                         signed_size_t const dim) {
    
    if (superKeyCompare(c->tuple, a->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_4_ab(KdNode<K,V,N>* const a,
                                      KdNode<K,V,N>* const b,
                                      KdNode<K,V,N>* const c,
                                      KdNode<K,V,N>* const d,
                                      signed_size_t const p,
                                      signed_size_t const dim) {
    
    if (superKeyCompare(c->tuple, d->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_1_4(KdNode<K,V,N>* const a,
                                   KdNode<K,V,N>* const b,
                                   KdNode<K,V,N>* const c,
                                   KdNode<K,V,N>* const d,
                                   signed_size_t const p,
                                   signed_size_t const dim) {
    
    if (superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_2_5_ab_cd(KdNode<K,V,N>* const a,
                                         KdNode<K,V,N>* const b,
                                         KdNode<K,V,N>* const c,
                                         KdNode<K,V,N>* const d,
                                         KdNode<K,V,N>* const e,
                                         signed_size_t const p,
                                         signed_size_t const dim) {
    
    if (superKeyCompare(c->tuple, a->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_2_5_ab(KdNode<K,V,N>* const a,
                                      KdNode<K,V,N>* const b,
                                      KdNode<K,V,N>* const c,
                                      KdNode<K,V,N>* const d,
                                      KdNode<K,V,N>* const e,
                                      signed_size_t const p,
                                      signed_size_t const dim) {
    
    if (superKeyCompare(c->tuple, d->tuple, p, dim) < 0) {
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
  static KdNode<K,V,N>* select_2_5(KdNode<K,V,N>* const a,
                                   KdNode<K,V,N>* const b,
                                   KdNode<K,V,N>* const c,
                                   KdNode<K,V,N>* const d,
                                   KdNode<K,V,N>* const e,
                                   signed_size_t const p,
                                   signed_size_t const dim) {
    
    if (superKeyCompare(a->tuple, b->tuple, p, dim) < 0) {
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
  static signed_size_t partition(KdNode<K,V,N>** const a,
                                 signed_size_t const start,
                                 signed_size_t const n,
                                 signed_size_t const size,
                                 signed_size_t const k,
                                 KdNode<K,V,N>** const medians,
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
        for (j = i; j > start && superKeyCompare(a[j - 1]->tuple, tmp->tuple, p, dim) > 0; --j) {
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
    //
    // NOTE, however, that NO value of INDEX_CUTOFF appears to improve
    // the performance of two threads relative to that of one thread.
    // Hence, the cost of spawning a child thread appears to exceed any
    // improvement in the performance of finding the index of the median
    // of medians that may be achieved via two threads.
#ifdef DUAL_THREAD_INDEX
    if (twoThreads && n > INDEX_CUTOFF) {

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
        throw runtime_error"\n\ncaught exception for index future in partition\n");
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
    if (superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p, dim) < 0) {
      ++i;
    }
    else if (superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p, dim) > 0) {
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
    if (superKeyCompare(a[start + i]->tuple, medianOfMedians->tuple, p, dim) > 0) {
      break;
    }
  }
#else
  // Search upward from the beginning of the array a in order to minimize the use of
  // the superKeyCompare function at the expense of greater use of the swap function.
  for (signed_size_t j = 0; j < n - 1; ++j) {
    if (superKeyCompare(a[start + j]->tuple, medianOfMedians->tuple, p, dim) < 0) {
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
static KdNode<K,V,N>* buildKdTree(KdNode<K,V,N>** reference,
                                  KdNode<K,V,N>** temporary,
                                  vector<signed_size_t> const& permutation,
                                  signed_size_t start, signed_size_t end, signed_size_t size,
                                  signed_size_t dim, signed_size_t maximumSubmitDepth,
                                  signed_size_t depth) {

  KdNode<K,V,N>* node = nullptr;

  // The partition permutes as x, y, z, w... and specifies the most significant key.
  signed_size_t const p = permutation[depth];

  if (end == start) {

    // Only one KdNode was passed to this method, so store it at this level of the tree.
    node = reference[start];

  }
  else if (end == start + 1) {

    // Two KdNodes were passed to this method in unsorted order, so store the
    // start KdNode at this level of the tree and determine whether to store the
    // end KdNode as the < child or the > child.
    node = reference[start];
    if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) > 0) {
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
    if (superKeyCompare(reference[start]->tuple, reference[mid]->tuple, p, dim) < 0) {
      // reference[start] < reference[mid]
      if (superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p, dim) < 0) {
        // reference[start] < reference[mid] < reference[end]
        node = reference[mid];
        node->ltChild = reference[start];
        node->gtChild = reference[end];
      }
      else {
        // reference[start] < reference[mid]; reference[end] < reference[mid]
        if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) < 0) {
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
      if (superKeyCompare(reference[start]->tuple, reference[end]->tuple, p, dim) < 0) {
        // reference[mid] < reference[start] < reference[end]
        node = reference[start];
        node->ltChild = reference[mid];
        node->gtChild = reference[end];
      }
      else {
        // reference[mid] < reference[start]; reference[end] < reference[start]
        if (superKeyCompare(reference[mid]->tuple, reference[end]->tuple, p, dim) < 0) {
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
      // The recursive call to buildKdTree must be placed in a lambda
      // expression because buildKdTree is a template not a function.
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
static KdNode<K,V,N>* buildKdTreePresorted(KdNode<K,V,N>** reference,
                                           KdNode<K,V,N>** temporary,
                                           vector<signed_size_t> const& permutation,
                                           signed_size_t start, signed_size_t end, signed_size_t size,
                                           signed_size_t dim, signed_size_t maximumSubmitDepth) {

  KdNode<K,V,N>* node = nullptr;

  // It is assumed that the reference array has been pre-sorted using the x:y:z:w... super key.
  signed_size_t const depth = 0;

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
 * The verifyKdTree function walks the k-d tree and checks that the
 * children of a node are in the correct branch of that node.
 *
 * Calling parameters:
 *
 * permutation - the permutation vector
 * dim - the number of dimensions
 * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
 * depth - the depth in the k-d tree
 *
 * returns: a count of the number of kdNodes in the k-d tree
 */
private:
signed_size_t verifyKdTree(vector<signed_size_t> const& permutation,
                           signed_size_t const dim,
                           signed_size_t const maximumSubmitDepth,
                           signed_size_t const depth) const {

  signed_size_t count = 1;

  // The partition cycles as x, y, z, w...
  signed_size_t const p = permutation[depth];

  if (ltChild != nullptr) {
    if (ltChild->tuple[p] > tuple[p]) {
      throw runtime_error("\n\nchild is > node in verifyKdTree\n");
    }
    if (superKeyCompare(ltChild->tuple, tuple, p, dim) >= 0) {
      throw runtime_error("\n\nchild is >= node in verifyKdTree\n");
    }
  }
  if (gtChild != nullptr) {
    if (gtChild->tuple[p] < tuple[p]) {
      throw runtime_error("\n\nchild is < node in verifyKdTree\n");
    }
    if (superKeyCompare(gtChild->tuple, tuple, p, dim) <= 0) {
      throw runtime_error("\n\nchild is <= node in verifyKdTree\n");
    }
  }

  // Verify the < branch with a child thread at as many levels of the tree as possible.
  // Create the child thread as high in the tree as possible for greater utilization.

  // Is a child thread available to verify the < branch?
  if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

    // No, so verify the < branch with the current thread.
    if (ltChild != nullptr) {
      count += ltChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
    }

    // Then verify the > branch with the current thread.
    if (gtChild != nullptr) {
      count += gtChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
    }
  }
  else {

    // Yes, so verify the < branch with a child thread. Note that a lambda
    // is required because this verifyKdTree function is not static. The use
    // of ref may be unnecessary in view of the [&] lambda argument specification.

    future<signed_size_t> verifyFuture;
    if (ltChild != nullptr) {
      verifyFuture =
        async(launch::async, [&] {
                               return ltChild->verifyKdTree(ref(permutation),
                                                            dim,
                                                            maximumSubmitDepth,
                                                            depth + 1);
                             });
    }

    // And simultaneously verify the > branch with the current thread.
    signed_size_t gtCount = 0;
    if (gtChild != nullptr) {
      gtCount = gtChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
    }

    // Wait for the child thread to finish execution.
    signed_size_t ltCount = 0;
    if (ltChild != nullptr) {
      try {
        ltCount = verifyFuture.get();
      }
      catch (exception const& e) {
        throw runtime_error("\n\ncaught exception for verify future in verifyKdTree\n");
      }
    }

    // Sum the counts returned by the child and current threads.
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
 * numDimensions - the number of dimensions
 * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
 *
 * returns: a KdNode pointer to the root of the k-d tree
 */
public:
static KdNode<K,V,N>* createKdTree(vector<KdNode<K,V,N>*>& kdNodes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const maximumSubmitDepth) {

  struct timespec startTime, endTime;

  // Create a temporary vector for use in sorting the kdNodes vector and building the k-d tree.
  startTime = getTime();
  vector<KdNode<K,V,N>*> temporary(kdNodes.size());

  // Sort the kdNodes vector using multiple threads. Importantly,
  // for compatibility with the 'permutation' vector initialized below,
  // use the first dimension (0) as the leading key of the super key.
  startTime = getTime();
  mergeSortReferenceAscending(kdNodes.data(), temporary.data(), 0, kdNodes.size() - 1,
                              0, numDimensions, maximumSubmitDepth, 0);
  endTime = getTime();
  double sortTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  // Remove KdNodes with duplicate tuples via one pass through the kdNodes vector.
  startTime = getTime();
  signed_size_t const end = removeDuplicates(kdNodes.data(), 0, numDimensions, kdNodes.size());
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

  // It is unnecessary to compute the partition coordinate upon each recursive call of
  // the buildKdTree or verifyKdTree functions because that coordinate depends only on
  // the depth of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, numDimensions, kdNodes.size());

  // Build the k-d tree with multiple threads if possible.
  auto const root = buildKdTreePresorted(kdNodes.data(), temporary.data(), permutation, 0, end,
                                         kdNodes.size(), numDimensions, maximumSubmitDepth);
  endTime = getTime();
  double const kdTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  // Verify the k-d tree and report the number of kdNodes.
  startTime = getTime();
  signed_size_t const numberOfNodes = root->verifyKdTree(permutation, numDimensions, maximumSubmitDepth, 0);
  endTime = getTime();
  double const verifyTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
  endTime = getTime();
  cout << "Number of nodes = " << numberOfNodes << endl << endl;

  cout << "totalTime = " << fixed << setprecision(2) << (sortTime + removeTime + kdTime + verifyTime)
       << "  sortTime = " << sortTime << "  removeTime = " << removeTime
       << "  kdTime = " << kdTime << "  verifyTime = " << verifyTime << endl << endl;

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
bool insideBounds(vector<K> const& queryLower,
                  vector<K> const& queryUpper,
                  vector<bool> const& enable) const {
    
  bool inside = true;
  for (size_t i = 0; i < queryLower.size(); ++i) {
    if (enable[i] && (queryLower[i] > tuple[i] || queryUpper[i] < tuple[i])) {
      inside = false;
      break;
    }
  }
  return inside;
}

/*
 * The regionSearch function searches the k-d tree recursively to find the KdNodes that
 * lie within a hyper-rectangle defined by the query lower and upper bounds.
 *
 * Calling parameters:
 *
 * result - a list<KdNode<K,V,N>*> that is passed by reference and modified
 * queryLower - the query lower bound vector
 * queryUpper - the query upper bound vector
 * permutation - vector that specifies permutation of the partition coordinate
 * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
 * depth - the depth in the k-d tree
 * enable - a vector that specifies the dimensions on which to prune the region search
 */
private:
void regionSearch(list<KdNode<K,V,N>*>& result,
                  vector<K> const& queryLower,
                  vector<K> const& queryUpper,
                  vector<signed_size_t> const& permutation,
                  signed_size_t const maximumSubmitDepth,
                  signed_size_t const depth,
                  vector<bool> const&  enable) {

  // The partition cycles as x, y, z, w...
  signed_size_t const p = permutation[depth];

  // If the KdNode is within the query hyper-rectangle for each of the k dimensions,
  // add the KdNode to the list of KdNodes that lie inside the hyper-cube. The
  // following loop is equivalent to the IN_REGION pseudo-Algol code proposed
  // by Jon Bentley in his CACM article.
  if (insideBounds(queryLower, queryUpper, enable)) {
    result.push_front(this);
  }

  // Determine whether to search the < and > branches of the k-d tree. Although
  // the superKeyCompare function can produce a different result for the == case
  // than does comparing only the leading keys of the super-keys, that result
  // will avoid unnecessary searching of a sub-tree (at the expense of a more
  // precise super-key comparison) but the unnecessary search appears not to
  // change the outcome of this recursive regionSearch function.
  //
  // Note that if the partition dimension is not enabled, both branches are searched.
#ifdef NO_SUPER_KEY
  bool const searchLT = ltChild != nullptr && (tuple[p] >= queryLower[p] || !enable[p]);
  bool const searchGT = gtChild != nullptr && (tuple[p] <= queryUpper[p] || !enable[p]);;
#else
  bool const searchLT = ltChild != nullptr && (superKeyCompare(tuple, queryLower.data(), p, queryLower.size()) >= 0
                                               || !enable[p]);
  bool const searchGT = gtChild != nullptr && (superKeyCompare(tuple, queryUpper.data(), p, queryLower.size()) <= 0
                                               || !enable[p]);
#endif

  // Do both branches require searching and is a child thread available?
  if (searchLT && searchGT && maximumSubmitDepth >= 0 && depth <= maximumSubmitDepth) {

    // Yes, both branches of the tree require searching and a child thread is available,
    // so prepare to search the < branch with a child thread.
    future<void> searchFuture;

    // Search the < branch?
    if (searchLT) {
        
      // Yes, search the < branch asynchronously with a child thread.
      // A lamba is required because this regionSearch function is not
      // static. The use of std::ref may be unnecessary in view of the
      // [&] lambda argument specification.
      list<KdNode<K,V,N>*> ltResult;
      searchFuture = async(launch::async, [&] {
                                            ltChild->regionSearch(ref(ltResult),
                                                                  ref(queryLower),
                                                                  ref(queryUpper),
                                                                  ref(permutation),
                                                                  maximumSubmitDepth,
                                                                  depth + 1,
                                                                  ref(enable));
                                          });
      // Search the > branch?
      list<KdNode<K,V,N>*> gtResult;
      if (searchGT) {
          
        // Yes, search the > branch  with the master thread.
        gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
      }

      // Get the result of searching the < branch with the child thread.
      try {
        searchFuture.get();
      }
      catch (exception const& e) {
        throw runtime_error("\n\ncaught exception for search future in regionSearch\n");
      }

      // Append the results of searching the < and > branches to the result (if any) for this KdNode.
      result.splice(result.end(), ltResult);
      result.splice(result.end(), gtResult);

    } else {

      // No, don't search the < branch. Search the > branch?
      list<KdNode<K,V,N>*> gtResult;
      if (searchGT) {
          
        // Yes, search the > branch  with the master thread.
        gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
      }

      // Append the result of searching the > branch to the result (if any) for this KdNode.
      result.splice(result.end(), gtResult);
    }

  } else {
      
    // No, both branches do not require searching. Search the < branch with the master thread?
    if (searchLT) {
      list<KdNode<K,V,N>*> ltResult;
      ltChild->regionSearch(ltResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
      result.splice(result.end(), ltResult);
    }

    // Search the > branch with the master thread?
    if (searchGT) {
      list<KdNode<K,V,N>*> gtResult;
      gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
      result.splice(result.end(), gtResult);
    }

  }
}

/*
 * Create a permutation vector.
 *
 * Calling parameters:
 *
 * permutation - the permutation vector that is passed by reference and modified
 * numDimensions - the number of dimensions
 * numCoordinates - the number of points in the coordinates vector
 */
private:
static
void createPermutation(vector<signed_size_t>& permutation,
                       signed_size_t const numDimensions,
                       signed_size_t const numCoordinates) {
    
  // Determine the maximum depth of the k-d tree, which is log2(numCoordinates).
  signed_size_t size = numCoordinates;
  signed_size_t maxDepth = 1;
  while (size > 0) {
    ++maxDepth;
    size >>= 1;
  }

  // Because the partition coordinate permutes in the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
  // (for e.g. 4-dimensional data), the leading key of the super key will be 0 at the
  // first level of the nascent tree, consistent with having sorted the reference array
  // using 0 as the leading key of the super key.
  permutation.resize(maxDepth);
  for (size_t i = 0; i < permutation.size(); ++i) {
    permutation[i] = i % numDimensions;
  }
}

/*
 * The searchRegion function searches the k-d tree to find the KdNodes that
 * lie within a hyper-rectangle defined by the query lower and upper bounds.
 *
 * Calling parameters:
 *
 * result - a list<KdNode<K,V,N>*> that is passed by reference and modified
 * queryLower - the query lower bound vector that is passed by reference and modified
 * queryUpper - the query upper bound vector that is passed by reference and modified
 * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
 * size - the number of points in the coordinates vector
 *
 * return a list of KdNodes that lie within the query hyper-rectangle
 */
public:
void searchRegion(list<KdNode<K,V,N>*>& result,
                  vector<K>& queryLower,
                  vector<K>& queryUpper,
                  signed_size_t const maximumSubmitDepth,
                  signed_size_t const size) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the regionSearch function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, queryLower.size(), size);
    
  // Ensure that each query lower bound <= the corresponding query upper bound.
  for (size_t i = 0; i < queryLower.size(); ++i) {
    if (queryLower[i] > queryUpper[i]) {
      auto const tmp = queryLower[i];
      queryLower[i] = queryUpper[i];
      queryUpper[i] = tmp;
    }
  }

  // Search the tree over all dimensions to obtain the resulting list of KdNodes.
  vector<bool> enable(queryLower.size(), true);
  regionSearch(result, queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
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
void searchRegion(list<KdNode<K,V,N>*>& result,
                  vector<K>& queryLower,
                  vector<K>& queryUpper,
                  signed_size_t const maximumSubmitDepth,
                  signed_size_t const size,
                  vector<bool> const& enable) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the regionSearch function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, queryLower.size(), size);
    
  // Ensure that each query lower bound <= the corresponding query upper bound.
  for (size_t i = 0; i < queryLower.size(); ++i) {
    if (queryLower[i] > queryUpper[i]) {
      auto const tmp = queryLower[i];
      queryLower[i] = queryUpper[i];
      queryUpper[i] = tmp;
    }
  }

  // Search the tree over the enabled dimensions to obtain the resulting list of KdNodes.
  regionSearch(result, queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
}

/*
 * Walk the k-d tree recursively and append to a list each KdNode that lies inside
 * the hyper-rectangle defined by the query lower and upper bounds.
 *
 * Calling parameters:
 *
 * result - a list of KdNode pointers that is passed by reference and modified
 * queryLower - the query lower bound vector
 * queryUpper - the query upper bound vector
 * enable - a vector that specifies the dimensions on which to test for insidedness
 *          and prune the search
 *
 * return a list of pointers to KdNodes that lie within the query hyper-rectangle.
 */
private:
void regionBrute(list<KdNode<K,V,N>*>& result,
                 vector<K> const& queryLower,
                 vector<K> const& queryUpper,
                 vector<bool> const& enable) {

  // Append the KdNode to the list if it lies inside the query bounds.
  if (insideBounds(queryLower, queryUpper, enable)) {
    result.push_front(this);
  }

  // Visit the < sub-tree.
  if (ltChild != nullptr) {
    list<KdNode<K,V,N>*> ltResult;
    ltChild->regionBrute(ltResult, queryLower, queryUpper, enable);
    result.splice(result.end(), ltResult);
  }

  // Visit the > sub-tree.
  if (gtChild != nullptr) {
    list<KdNode<K,V,N>*> gtResult;
    gtChild->regionBrute(gtResult, queryLower, queryUpper, enable);
    result.splice(result.end(), gtResult);
  }
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
void bruteRegion(list<KdNode<K,V,N>*>& result,
                 vector<K>& queryLower,
                 vector<K>& queryUpper) {

  // Search over all dimensions.
  vector<bool> enable(queryLower.size(), true);

  // Ensure that each query lower bound <= the corresponding query upper bound.
  for (size_t i = 0; i < queryLower.size(); ++i) {
    if (queryLower[i] > queryUpper[i]) {
      auto const tmp = queryLower[i];
      queryLower[i] = queryUpper[i];
      queryUpper[i] = tmp;
    }
  }

  // Walk the k-d tree recursively.
  regionBrute(result, queryLower, queryUpper, enable);
}

/*
 * Search the k-d tree recursively for up to M nearest geometric neighbors to a query point
 * by adding them to the NearestNeighborHeap that stores up to M neighbors.  Exclude from the
 * search any branch of the tree wherein it is guaranteed that all nodes in that branch are
 * farther away than the current farthest node stored in the heap. Details of the search
 * algorithm are described by Friedman et al. in "An Algorithm for Finding Best Matches in
 * Logarithmic Expected Time", ACM Transactions on Mathematical Software, 3: 209-226, 1977.
 *
 * Calling parameters:
 *
 * heap - an instance of NearestNeighborHeap that is built relative to a query point
 * permutation - vector that specifies permutation of the partition coordinate
 * depth - depth in the k-d tree
 */
private:
void nearestNeighbors(NearestNeighborHeap<K,V,N>& heap,
                      vector<signed_size_t> const& permutation,
                      signed_size_t const depth) {

  // The partition permutes as x, y, z, w...
  signed_size_t p = permutation[depth];

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
    double const dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
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
    double const dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
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
void findNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighbors,
                          vector<K> const& query,
                          signed_size_t const numNeighbors,
                          signed_size_t const size) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the nearestNeighbors function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, query.size(), size);

  // Create the heap and search the k-d tree for nearest neighbors.
  NearestNeighborHeap<K,V,N> heap(query, numNeighbors);
  nearestNeighbors(heap, permutation, 0);

  // Empty the heap by successively removing the top of the heap and prepending it to a list.
  // Remove only the number of heap entries present.
  signed_size_t const heapDepth = heap.heapDepth();
  for (signed_size_t i = 0; i < heapDepth; ++i) {
    neighbors.push_front(heap.removeTop());
  }
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
void findNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighbors,
                          vector<K> const& query,
                          signed_size_t const numNeighbors,
                          signed_size_t const size,
                          vector<bool> const& enable) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the nearestNeighbors function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, query.size(), size);

  // Create the heap and search the k-d tree for nearest neighbors.
  NearestNeighborHeap<K,V,N> heap(query, numNeighbors, enable);
  nearestNeighbors(heap, permutation, 0);

  // Empty the heap by successively removing the top of the heap and prepending it to a list.
  // Remove only the number of heap entries present.
  signed_size_t const heapDepth = heap.heapDepth();;
  for (signed_size_t i = 0; i < heapDepth; ++i) {
    neighbors.push_front(heap.removeTop());
  }
}
  
/*
 * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
 *
 * Calling parameters:
 *
 * neighbors - the nearest neighbors list that is passed by reference and modified.
 * query - the query vector
 * permutation - vector that specifies permutation of the partition coordinate
 * numNeighbors - the number M of nearest neighbors to attempt to find
 */
public:
void findNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighbors,
                          vector<K> const& query,
                          vector<signed_size_t> const& permutation,
                          signed_size_t const numNeighbors) {
    
  // Create the heap and search the k-d tree for nearest neighbors.
  NearestNeighborHeap<K,V,N> heap(query, numNeighbors);
  nearestNeighbors(heap, permutation, 0);

  // Empty the heap by successively removing the top of the heap and prepending it to a list.
  // Remove only the number of heap entries present.
  signed_size_t const heapDepth = heap.heapDepth();;
  for (signed_size_t i = 0; i < heapDepth; ++i) {
    neighbors.push_front(heap.removeTop());
  }
}

/*
 * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
 *
 * Calling parameters:
 *
 * neighbors - the nearest neighbors list that is passed by reference and modified.
 * query - the query vector
 * permutation - vector that specifies permutation of the partition coordinate
 * numNeighbors - the number M of nearest neighbors to attempt to  find
 * enable - a vector that specifies the dimensions for which to test distance
 */
public:
void findNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighbors,
                          vector<K> const& query,
                          vector<signed_size_t> const& permutation,
                          signed_size_t const numNeighbors,
                          vector<bool> const& enable) {
    
  // Create the heap and search the k-d tree for nearest neighbors.
  NearestNeighborHeap<K,V,N> heap(query, numNeighbors, enable);
  nearestNeighbors(heap, permutation, 0);

  // Empty the heap by successively removing the top of the heap and prepending it to a list.
  // Remove only the number of heap entries present.
  signed_size_t const heapDepth = heap.heapDepth();;
  for (signed_size_t i = 0; i < heapDepth; ++i) {
    neighbors.push_front(heap.removeTop());
  }
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
void verifyNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighborsFast,
                            forward_list< pair<double, KdNode<K,V,N>*> >& neighborsSlow) const {
    
  auto itf1 = neighborsFast.begin();
  auto its1 = neighborsSlow.begin();
    
  // Compare the first k-d tree (fast) distance to the first brute-force (slow) distance.
  if (itf1->first != its1->first) {
    ostringstream buffer;
    buffer << "\n\nfast distance[0] = " << itf1->first << "  !=  slow distance[0] = " << its1->first << endl;
    throw runtime_error(buffer.str());
  }
    
  // Compare the first k-d tree KdNode pointer to the first brute-force KdNode pointer.
  if (itf1->second != its1->second) {
    throw runtime_error("\n\nfast KdNode*[0] != slow KdNode*[0]\n");
  }

  // Compare the remaining distances and KdNode*
  auto itf2 = itf1;
  auto its2 = its1;
  ++itf2;
  ++its2;
  signed_size_t i = 1;
  for ( ; itf2 != neighborsFast.end(); ++itf1, ++its1, ++itf2, ++its2, ++i) {
    // Ensure that the fast distances increase monotonically.
    if (itf1->first > itf2->first) {
      ostringstream buffer;
      buffer << "\n\nfast distance[" << (i-1) << "] = " << itf1->first << "  >  fast distance[" << i << "] = " << itf2->first << endl;
      throw runtime_error(buffer.str());
    }
    // Ensure that the slow distances increase monotonically.
    if (its1->first > its2->first) {
      ostringstream buffer;
      buffer << "\n\nslow distance[" << (i-1) << "] = " << its1->first << "  >  slow distance[" << i << "] = " << its2->first << endl;
      throw runtime_error(buffer.str());
    }
    // Compare the ith k-d tree distance to the ith brute-force distance.
    if (itf2->first != its2->first) {
      ostringstream buffer;
      buffer << "\n\nfast distance[" << i << "] = " << itf2->first << "  !=  slow distance[" << i << "] = " << its2->first << endl;
      throw runtime_error(buffer.str());
    }
    // Compare the ith k-d tree KdNode pointer to the ith brute-force KdNode pointer.
    if (itf2->second != its2->second) {
      ostringstream buffer;
      buffer << "\n\nfast KdNode*[" << i << "]  !=  slow KdNode*[" << i << "]\n";
      throw runtime_error(buffer.str());
    }
  }
}
  
/*
 * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
 * and add those nearest neighbors to a reverse nearest neighbors map and to a
 * nearest neighbors map.
 *
 * Calling parameters:
 *
 * nn - the nearest neighbors map that is passed by reference and modified
 * rnn - the reverse nearest neighbors map that is passed by reference and modified
 * permutation - vector that specifies permutation of the partition coordinate
 * root - the root of the k-d tree where a search for nearest neighbors must begin
 * numDimensions - the dimensionality k of the k-d tree
 * numNeighbors - the number M of nearest neighbors to attempt to find
 */
private:
void reverseNearestNeighbors(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& nn,
                             map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
                             vector<signed_size_t> const& permutation,
                             KdNode<K,V,N>* const root,
                             signed_size_t const numDimensions,
                             signed_size_t const numNeighbors) {

  // Create a query point from the KdNode's tuple, find at most the M
  // nearest neighbors to it, prepend those neighbors to a nearest
  // neighbors list, remove the first element of the list (which is
  // the query KdNode), and store a pointer to the list in the nn map.
  vector<K> const query(tuple, tuple + numDimensions);
  auto nnList = nn[this];
  root->findNearestNeighbors(*nnList, query, permutation, numNeighbors);
  nnList->pop_front();

  // Iterate over the remaining list of nearest neighbors and prepend
  // the query KdNode to the reverse nearest neighbors list at the
  // map entry for the nearest neighbor. There is no need to update
  // the map because it contains a pointer to a list, not a list.
  for (auto it = nnList->begin(); it != nnList->end(); ++it) {
    rnn[it->second]->push_front(make_pair(it->first, this));
  }
    
  // Visit the < sub-tree.
  if (ltChild != nullptr) {
    ltChild->reverseNearestNeighbors(nn, rnn, permutation, root, numDimensions, numNeighbors);
  }
    
  // Visit the > sub-tree.
  if (gtChild != nullptr) {
    gtChild->reverseNearestNeighbors(nn, rnn, permutation, root, numDimensions, numNeighbors);
  }
}

/*
 * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
 * and add those nearest neighbors to a reverse nearest neighbors map and to a
 * nearest neighbors map.
 *
 * Calling parameters:
 *
 * nn - the nearest neighbors map that is passed by reference and modified
 * rnn - the reverse nearest neighbors map that is passed by reference and modified
 * permutation - vector that specifies permutation of the partition coordinate
 * root - the root of the k-d tree where a search for nearest neighbors must begin
 * numDimensions - the dimensionality k of the k-d tree
 * numNeighbors - the number M of nearest neighbors to attempt to find
 * enable - a vector that specifies the dimensions for which to test distance
 */
private:
void reverseNearestNeighbors(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& nn,
                             map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
                             vector<signed_size_t> const& permutation,
                             KdNode<K,V,N>* const root,
                             signed_size_t const numDimensions,
                             signed_size_t const numNeighbors,
                             vector<bool> const& enable) {

  // Create a query point from the KdNode's tuple, find at most the M
  // nearest neighbors to it, prepend those neighbors to a nearest
  // neighbors list, remove the first element of the list (which is
  // the query KdNode), and store a pointer to the list in the nn map.
  vector<K> const query(tuple, tuple + numDimensions);
  auto nnList = nn[this];
  root->findNearestNeighbors(*nnList, query, permutation, numNeighbors, enable);
  nnList->pop_front();

  // Iterate over the remaining list of nearest neighbors and prepend
  // the query KdNode to the reverse nearest neighbors list at the
  // map entry for the nearest neighbor. There is no need to update
  // the map because it contains a pointer to a list, not a list.
  for (auto it = nnList->begin(); it != nnList->end(); ++it) {
    rnn[it->second]->push_front(make_pair(it->first, this));
  }
    
  // Visit the < sub-tree.
  if (ltChild != nullptr) {
    ltChild->reverseNearestNeighbors(nn, rnn, permutation, root, numDimensions, numNeighbors, enable);
  }
    
  // Visit the > sub-tree.
  if (gtChild != nullptr) {
    gtChild->reverseNearestNeighbors(nn, rnn, permutation, root, numDimensions, numNeighbors, enable);
  }
}

/*
 * Walk the k-d tree and initialize the reverse nearest neighbors map.
 *
 * Calling parameters:
 *
 * rnn - the reverse nearest neighbors map that is called by reference and modified
 * rnnLists - a vector of reverse nearest neighbors lists
 * index - an index into the vector that is passed by reference and modified
 */
private:
void initMap(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
             vector< forward_list< pair<double, KdNode<K,V,N>*> > >& rnnLists,
             size_t& index) {

  // Create a map entry that points to the vector element (which is a list)
  // and increment the vector index to the next vector element.
  rnn[this] = &rnnLists[index];
  ++index;
    
  // Visit the < sub-tree.
  if (ltChild != nullptr) {
    ltChild->initMap(rnn, rnnLists, index);
  }
    
  // Visit the > sub-tree.
  if (gtChild != nullptr) {
    gtChild->initMap(rnn, rnnLists, index);
  }
}
  
/*
 * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
 * and add those nearest neighbors to a reverse nearest neighbors map and to a
 * nearest neighbors map.
 *
 * Each element of the nearest neighbors map contains a query KdNode and a list
 * of KdNodes that are nearest neighbors to the query KdNode.
 *
 * Each element of the reverse nearest neighbors map contains a reference KdNode
 * and a list of KdNodes to which the reference KdNode is a nearest neighbor.
 * The concept of reverse nearest neighbors was first described by F. Korn and
 * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
 * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
 * Data, pages 201-212.
 *
 * Calling parameters:
 *
 * nn - the nearest neighbors map that is passed by reference and modified
 * rnn - the reverse nearest neighbors map that is passed by reference and modified
 * nnLists - a vector of nearest neighbors lists
 * rnnLists - a vector of reverse nearest neighbors lists
 * numDimensions - the dimensionality k of the k-d tree
 * numNeighbors - the number M of nearest neighbors to attempt to find
 */
public:
void findReverseNearestNeighbors(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& nn,
                                 map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
                                 vector< forward_list< pair<double, KdNode<K,V,N>*> > >& nnLists,
                                 vector< forward_list< pair<double, KdNode<K,V,N>*> > >& rnnLists,
                                 signed_size_t const numDimensions,
                                 signed_size_t const numNeighbors) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the nearestNeighbors function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, numDimensions, rnnLists.size());

  // Initialize the nearest neighbors and reverse nearest neighbors maps.
  size_t index = 0;
  initMap(nn, nnLists, index);
  index = 0;
  initMap(rnn, rnnLists, index);
    
  // Walk the k-d tree and build the reverse nearest neighbors lists.
  reverseNearestNeighbors(nn, rnn, permutation, this, numDimensions, numNeighbors);
}

/*
 * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
 * and add those nearest neighbors to a reverse nearest neighbors map and to a
 * nearest neighbors map.
 *
 * Each element of the nearest neighbors map contains a query KdNode and a list
 * of KdNodes that are nearest neighbors to the query KdNode.
 *
 * Each element of the reverse nearest neighbors map contains a reference KdNode
 * and a list of KdNodes to which the reference KdNode is a nearest neighbor.
 * The concept of reverse nearest neighbors was first described by F. Korn and
 * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
 * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
 * Data, pages 201-212.
 *
 * Calling parameters:
 *
 * nn - the nearest neighbors map that is passed by reference and modified
 * rnn - the reverse nearest neighbors map that is passed by reference and modified
 * nnLists - a vector of nearest neighbors lists
 * rnnLists - a vector of reverse nearest neighbors lists
 * numDimensions - the dimensionality k of the k-d tree
 * numNeighbors - the number M of nearest neighbors to attempt to find
 * enable - a vector that specifies the dimensions for which to test distance
 */
public:
void findReverseNearestNeighbors(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& nn,
                                 map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
                                 vector< forward_list< pair<double, KdNode<K,V,N>*> > >& nnLists,
                                 vector< forward_list< pair<double, KdNode<K,V,N>*> > >& rnnLists,
                                 signed_size_t const numDimensions,
                                 signed_size_t const numNeighbors,
                                 vector<bool> const& enable) {
    
  // It is unnecessary to compute the partition coordinate upon each recursive call
  // of the nearestNeighbors function because that coordinate depends only on the depth
  // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
  vector<signed_size_t> permutation;
  createPermutation(permutation, numDimensions, rnnLists.size());

  // Initialize the nearest neighbors and reverse nearest neighbors maps.
  size_t index = 0;
  initMap(nn, nnLists, index);
  index = 0;
  initMap(rnn, rnnLists, index);
    
  // Walk the k-d tree and build the reverse nearest neighbors lists.
  reverseNearestNeighbors(nn, rnn, permutation, this, numDimensions, numNeighbors, enable);
}

/*
 * Verify the correctness of the reverse nearest neighbors map.
 *
 * Calling parameter:
 *
 * nn - the nearest neighbors map
 * rnn - the reverse nearest neighbors map.
 *
 * Although this function does not directly access the k-d tree, it requires the persistence
 * of the k-d tree for access to the KdNodes via the maps. Hence, this function is not static.
 */
void verifyReverseNeighbors(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& nn,
                            map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn) const {

  // Iterate through the reverse nearest neighbors map and if a list is not empty,
  // verify the correctness of each list.
  for (auto it = rnn.begin(); it != rnn.end(); ++it) {
    auto const rnnListPtr = it->second;
    if (!rnnListPtr->empty()) {
      // Get the KdNode that is a nearest neighbor to all KdNodes on the list and
      // verify that it is indeed a nearest neighbor to each KdNode on the list.
      auto const kdNodePtr = it->first;
      for (auto listIt = rnnListPtr->begin(); listIt != rnnListPtr->end(); ++listIt) {
        auto const kdPtr = listIt->second;
        // Get the nearest neighbor list for the KdNode from the nearest neighbors map
        // and verify that the list contains KdNodePtr.
        auto nnListPtr = nn[kdPtr];
        bool match = false;
        for (auto nnIt = nnListPtr->begin(); nnIt != nnListPtr->end(); ++nnIt) {
          if (kdNodePtr == nnIt->second) {
            match = true;
            break;
          }
        }
        if (!match) {
          throw runtime_error("\n\nnode is not a nearest neighbor\n");
        }
      }
    }
  }
}

/*
 * Calculate the mean and standard deviation of the distances and list sizes in a map.
 *
 * Calling parameter:
 *
 * nmap - a map from KdNode pointer to list.
 *
 * Although this function does not directly access the k-d tree, it requires the persistence
 * of the k-d tree for access to the KdNodes via the map. Hence, this function is not static.
 */
void calculateMeanStd(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn,
                      double& meanSize, double& stdSize, double& meanDist, double& stdDist) const {

  // Count the number of map entries that have non-empty lists
  // and sum the distances and list lengths.
  size_t count = 0;
  double sumDist = 0.0, sumDist2 = 0.0, sumSize = 0.0, sumSize2 = 0.0;
  for (auto mapIt = rnn.begin(); mapIt != rnn.end(); ++mapIt) {
    auto const rnnListPtr = mapIt->second;
    if (!rnnListPtr->empty()) {
      ++count;
      double const size = static_cast<double>(distance(rnnListPtr->begin(), rnnListPtr->end()));
      sumSize += size;
      sumSize2 += size * size;
      for (auto listIt = rnnListPtr->begin(); listIt != rnnListPtr->end(); ++listIt) {
        double const dist2 = listIt->first;
        sumDist += sqrt(dist2);
        sumDist2 += dist2;
      }
    }
  }
  double mapSize = static_cast<double>(count);
  meanSize = sumSize / mapSize;
  stdSize = sqrt((sumSize2 - (sumSize * sumSize / mapSize)) / (mapSize - 1.0));
  meanDist = sumDist / sumSize;
  stdDist = sqrt((sumDist2 - (sumDist * sumDist / sumSize)) / (sumSize - 1.0));
}

/*
 * Count the number of non-empty lists in a map.
 *
 * Calling parameter:
 *
 * rnn - the map
 *
 * Although this function does not directly access the k-d tree, it requires the persistence
 * of the k-d tree for access to the KdNodes via the map. Hence, this function is not static.
 */
public:
size_t nonEmptyLists(map< KdNode<K,V,N>*, forward_list< pair<double, KdNode<K,V,N>*> >* >& rnn) const {

  size_t count = 0;
  for (auto it = rnn.begin(); it != rnn.end(); ++it) {
    if(!it->second->empty()) {
      ++count;
    }
  }
  return count;
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
 * neighbors - the nearest neighbors list that is passed by reference and modified.
 * query - the query vector
 * numNeighbors - the number M of nearest neighbors to find
 */
public:
void bruteNearestNeighbors(forward_list< pair<double, KdNode<K,V,N>*> >& neighbors,
                           vector<K> const& query,
                           signed_size_t const numNeighbors) {
    
  // Create the heap, walk the k-d tree, and attempt to add each KdNode to the heap.
  NearestNeighborHeap<K,V,N> heap(query, numNeighbors);
  allNeighbors(heap);

  // Empty the heap by successively removing the top of the heap and appending it to a list.
  for (signed_size_t i = 0; i < numNeighbors; ++i) {
    neighbors.push_front(heap.removeTop());
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
    
  cout << "(" << tuple[0] << ",";
  for (signed_size_t i = 1; i < dim - 1; ++i) cout << tuple[i] << ",";
  cout << tuple[dim - 1] << ")";
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
    
  cout << "(" << tuple[0] << ",";
  for (size_t i = 1; i < tuple.size() - 1; ++i) cout << tuple[i] << ",";
  cout << tuple[tuple.size() - 1] << ")";
}

/*
 * The printTuples func prints all tuples in a list.
 *
 * Calling parameters:
 *
 * regionList - a list of KdNodes returned by a region search
 * maximumNumberOfNodesToPrint - the maximum number of KdNodes to print
 * numDimensions - the number of dimensions
 *
 * Because this function does not access the k-d tree, it could be static.
 * However, calling it as a static function requires speicification of a
 * type, so calling it as a non-static function is less cumbersome.
 */
public:
void printTuples(list< KdNode<K,V,N>* > const& regionList,
                 signed_size_t const maximumNumberOfNodesToPrint,
                 signed_size_t const numDimensions) const {
    
  if (regionList.size() != 0) {
    signed_size_t maxNodesToPrint = maximumNumberOfNodesToPrint;
    for (auto it = regionList.begin(); it != regionList.end(); ++it) {
      printTuple((*it)->getTuple(), numDimensions);
      cout << endl;
      --maxNodesToPrint;
      if (maxNodesToPrint == 0) {
        break;
      }
    }
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
void printKdTree(signed_size_t const dim,
                 signed_size_t const depth) const {
  if (gtChild != nullptr) {
    gtChild->printKdTree(dim, depth + 1);
  }
  for (signed_size_t i = 0; i < depth; ++i) cout << "       ";
  printTuple(tuple, dim);
  cout << endl;
  if (ltChild != nullptr) {
    ltChild->printKdTree(dim, depth + 1);
  }
}
}; // class KdNode

/*
 * The NearestNeighborHeap Class implements a fixed length heap of both containing both a KdNode and Euclidean distance
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
  vector<K> query; // query point for which nearest neighbors will be found
  vector<bool> enable;
private:
  signed_size_t reqDepth; // requested number of nearest neighbors
  vector<KdNode<K,V,N>* > nodes; // vector of pointers to KdNodes that are the nearest neighbors
  vector<double> dists; // vector of squared distances
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
  NearestNeighborHeap(vector<K> const& query, signed_size_t numNeighbors) {
    this->nodes.resize(numNeighbors + 1, nullptr); // heap of KdNode* (address 0 is unused)
    this->dists.resize(numNeighbors + 1, 0); // corresponding heap of distances (initialized to 0)
    this->reqDepth = numNeighbors;
    this->curDepth = 0;
    this->query = query;
    this->enable.assign(query.size(), true);
  }
  
  /*
   * Constructor that enables distance test for only specified dimensions
   *
   * Calling parameters:
   *
   * query - a vector that defines the query point
   * numNeighbors - the number of nearest neighbors desired
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  NearestNeighborHeap(vector<K> const& query, signed_size_t numNeighbors, vector<bool> const& enable) {
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
  void swap(signed_size_t const i,
            signed_size_t const j) {
    
    double const tempDist = dists[i];
    auto const tempNode = nodes[i];
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
   * kk - the index of the element
   */
private:
  void rise(signed_size_t const kk) {

    signed_size_t k = kk;
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
   * kk - the index of the element
   */
private:
  void fall(signed_size_t const kk) {

    signed_size_t k = kk;
    while (2*k <= curDepth) {
      signed_size_t j = 2*k;
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
  pair<double, KdNode<K,V,N>*> removeTop() {
    pair<double, KdNode<K,V,N>*> returnPair = make_pair(dists[1], nodes[1]);
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
  void add(KdNode<K,V,N>* const node) {
    // Find the distance by subtracting the query from the tuple and
    // calculating the sum of the squared distances. Note that conversion
    // from type T to double may result in loss of precision but avoids
    // the possibility of integer overflow.
    double dist2 = 0.0;
    for (size_t i = 0; i < query.size(); ++i) {
      // Add the squared coordinate distance only if the dimension is enabled.
      if (enable[i]) {
        K const  comp = node->tuple[i] - query[i];
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

  /* Return the current depth of the heap, i.e., the number of nearest nodes/distances elements on the heap. */
public:
  signed_size_t heapDepth() {
    return curDepth;
  }
}; // class NearestNeighborHeap

/*
 * The randomLongInInterval function creates a random kdKey_t in the interval [min, max].  See
 * http://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
 *
 * Calling parameters:
 *
 * min - the minimum kdKey_t value desired
 * max - the maximum kdKey_t value desired
 *
 * returns: a random kdKey_t
 */
static kdKey_t randomLongInInterval(kdKey_t const min,
                                    kdKey_t const max) {
  
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
  signed_size_t maximumNumberOfNodesToPrint = 5;
  kdKey_t searchDistance = 1000000000000000000L;
  bool bruteForceSearch = false;
  bool bruteForceRegion = false;
  bool reverseNearestNeighbors = false;

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
    if (0 == strcmp(argv[i], "-r") || 0 == strcmp(argv[i], "--reverseNearestNeighbors")) {
      reverseNearestNeighbors = !reverseNearestNeighbors;
      continue;
    }
    {
      ostringstream buffer;
      buffer << "\n\nillegal command-line argument: " << argv[i] << endl;
      throw runtime_error(buffer.str());
    }
  }

  // Declare and initialize the kdNodes vector and initialize it with tuples,
  // for example (x,y,z,w), in the half-open interval [0, MAX] where MAX is the
  // maximum value for the kdKey_t type. Create extraPoints-1 duplicate elements,
  // where extraPoints <= numPoints, to test the removal of duplicate points.
  extraPoints = (extraPoints <= numPoints) ? extraPoints : numPoints;
  std::vector<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> kdNodes(numPoints + extraPoints - 1);
  // Create each KdNode with a tuple (all elements == 0) and a string value
  // obtained from the index into the kdNodes vector. The string value is
  // a proxy for a KdNode label.
  for (size_t i = 0; i < kdNodes.size(); ++i) {
    ostringstream buffer;
    buffer << i;
    kdNodes[i] = new KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>(buffer.str());
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
  cout << endl << "Max number of threads = " << numThreads << "  max submit depth = "
       << maximumSubmitDepth << endl << endl;

  // Create the k-d tree.
  auto const root = KdNode<kdKey_t, kdValue_t,K_DIMENSIONALITY>::createKdTree(kdNodes, K_DIMENSIONALITY, maximumSubmitDepth);

  // Search the k-d tree via region search for the KdNodes that lie within a hyper-cube centered near the origin.
  vector<kdKey_t> query(K_DIMENSIONALITY);
  vector<kdKey_t> queryLower(K_DIMENSIONALITY);
  vector<kdKey_t> queryUpper(K_DIMENSIONALITY);
  for (signed_size_t i = 0; i < K_DIMENSIONALITY; ++i) {
    query[i] = i;
    queryLower[i] = query[i] - searchDistance;
    queryUpper[i] = query[i] + searchDistance;
  }
  startTime = getTime();
  list<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> regionFast;
  root->searchRegion(regionFast, queryLower, queryUpper, maximumSubmitDepth, kdNodes.size());
  endTime = getTime();
  double const fastRegionTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  cout << "fast region time = " << fixed << setprecision(6) << fastRegionTime << " seconds" << endl << endl;

  cout << regionFast.size() << " nodes within " << searchDistance << " units of ";
  root->printTuple(query);
  cout << " in all dimensions." << endl << endl;
  cout << "List of the first <= " << maximumNumberOfNodesToPrint << " fast search k-d nodes within a "
       << searchDistance << "-unit search distance follows:" << endl << endl;
  regionFast.sort();
  root->printTuples(regionFast, maximumNumberOfNodesToPrint, K_DIMENSIONALITY);
  cout << endl;

  // Verify that no duplicate KdNodes exist on the list returned from region search.
  auto itr1 = regionFast.begin();
  auto itr2 = itr1;
  ++itr2;
  for ( ; itr2 != regionFast.end(); ++itr1, ++itr2) {
    if (*itr1 == *itr2) {
      throw runtime_error("\n\nduplicate KdNode* on region-search list\n");
    }
  }

  // Search the k-d tree via brute force for the KdNodes that lie within a hyper-cube centered near the origin.
  if (bruteForceRegion) {
    startTime = getTime();
    list<KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> regionSlow;
    root->bruteRegion(regionSlow, queryLower, queryUpper);
    endTime = getTime();
    double const slowRegionTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    cout << "slow region time = " << fixed << setprecision(6) << slowRegionTime << " seconds" << endl << endl;

    cout << regionSlow.size() << " nodes within " << searchDistance << " units of ";
    root->printTuple(query);
    cout << " in all dimensions." << endl << endl;
    cout << "List of the first <= " << maximumNumberOfNodesToPrint << " slow search k-d nodes within a "
         << searchDistance << "-unit search distance follows:" << endl << endl;
    regionSlow.sort();
    root->printTuples(regionSlow, maximumNumberOfNodesToPrint, K_DIMENSIONALITY);
    cout << endl;

    // Verify that the region-search and brute-force lists are identical.
    // Both lists must be sorted before the KdNode* comparisons are
    // performed below because the region search and brute-force search
    // algorithms do not produce lists wherein KdNodes are prepended
    // in the same order.
    auto itrf = regionFast.begin();
    for (auto itrs = regionSlow.begin(); itrs != regionSlow.end(); ++itrf, ++itrs) {
      if (*itrf != *itrs) {
        throw runtime_error("\n\nnon-identical region-search and brute-force lists\n");
      }
    }
  }

  // It is impossible to find more nearest neighbors than there are points.
  numNeighbors = min(numNeighbors, numPoints + extraPoints + 1);

  // Search the k-d tree for up to numNeighbors nearest neighbors to the first tuple.
  startTime = getTime();
  forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> > neighborsFast;
  root->findNearestNeighbors(neighborsFast, query, numNeighbors, kdNodes.size());
  endTime = getTime();
  double const fastNeighborTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  cout << "fast neighbor time = " << fixed << setprecision(6) << fastNeighborTime << " seconds" << endl << endl;
  cout << "fast neighbor list size = " << distance(neighborsFast.begin(), neighborsFast.end()) << endl << endl;

  // Find nearest neighbors via brute force if requested.
  if (bruteForceSearch) {
    startTime = getTime();
    forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> > neighborsSlow;
    // Find only the number of nearest neighbors returned by findNearestNeighbors above.
    root->bruteNearestNeighbors(neighborsSlow, query, distance(neighborsFast.begin(), neighborsFast.end()));
    endTime = getTime();
    double const slowNeighborTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    cout << "slow neighbor time = " << fixed << setprecision(6) << slowNeighborTime << " seconds" << endl << endl;
    cout << "slow neighbor list size = " << distance(neighborsSlow.begin(), neighborsSlow.end()) << endl << endl;

    // Verify the consistency between the nearest neighbors lists
    // found by k-d tree search and by brute force.
    root->verifyNearestNeighbors(neighborsFast, neighborsSlow);

    // Print the fast and slow distances squared.
    cout << "fast and slow distances squared follow:" << endl << endl;
    auto itf = neighborsFast.begin();
    auto its = neighborsSlow.begin();
    for ( ; itf != neighborsFast.end(); ++itf, ++its) {
      cout << fixed << setprecision(0) << itf->first << "\t" << its->first << endl;
    }
    cout << endl;
  }

  // Optionally construct a nearest neighbor map and a reverse nearest neighbors map.
  // Each map element contains a list. The lists are allocated within vectors below
  // to permit the lists to be manipulated within the maps via pointers to lists, which
  // avoids the possibility of copying lists that may be time-consuming.
  if (reverseNearestNeighbors) {
    startTime = getTime();
    map< KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*,
         forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> >* > nn;
    map< KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*,
         forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> >* > rnn;
    vector< forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> > > nnLists(kdNodes.size());
    vector< forward_list< pair<double, KdNode<kdKey_t, kdValue_t, K_DIMENSIONALITY>*> > > rnnLists(kdNodes.size());
    root->findReverseNearestNeighbors(nn, rnn, nnLists, rnnLists, K_DIMENSIONALITY, numNeighbors);
    endTime = getTime();
    double const reverseNearestNeighborTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    cout << "reverse nearest neighbor time = " << fixed << setprecision(2) << reverseNearestNeighborTime << " seconds" << endl << endl;
    cout << "number of non-empty nearest neighbors lists = " << root->nonEmptyLists(nn) << endl;
    cout << "number of non-empty reverse nearest neighbors lists = " << root->nonEmptyLists(rnn) << endl << endl;

    // Report the mean and standard deviation distance and number of reverse nearest neighbors.
    double meanSize, stdSize, meanDist, stdDist;
    root->calculateMeanStd(rnn, meanSize, stdSize, meanDist, stdDist);
    cout << "mean reverse nearest neighbor distance = " << scientific << meanDist
         << "  standard deviation = " << stdDist << endl;
    cout << "mean reverse nearest neighbor list size = " << fixed << setprecision(3) << meanSize
         << "  standard deviation = " << stdSize << endl << endl;

    // Verify the consistency between the nearest neighbors and reverse nearest neighbors maps.
    startTime = getTime();
    root->verifyReverseNeighbors(nn, rnn);
    endTime = getTime();
    double verifyReverseTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

    cout << "verify reverse nearest neighbor time = " << fixed << setprecision(2) << verifyReverseTime << " seconds" << endl << endl;
  }

  // Delete the k-d tree.
  startTime = getTime();
  root->deleteKdTree();
  endTime = getTime();
  double const deleteTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

  cout << "deleteTime = " << fixed << setprecision(6) << deleteTime << " seconds" << endl << endl;

  return 0;
}
#endif // #ifdef TEST_KD_TREE
