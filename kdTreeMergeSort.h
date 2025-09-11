/*
 * Copyright (c) 2015, 2021, 2023, 2024 Russell A. Brown
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

#ifndef KD_TREE_MERGE_SORT_H
#define KD_TREE_MERGE_SORT_H

/* A cutoff for switching from merge sort to insertion sort in KdNode::mergeSort* */
#ifndef INSERTION_SORT_CUTOFF
#define INSERTION_SORT_CUTOFF 15
#endif

/* A cutoff for switching from 1 to 2 threads to merge reference arrays in KdNode::mergeSort* */
#ifndef MERGE_CUTOFF
#define MERGE_CUTOFF 4096
#endif

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

/* The merge sort functions */
template <typename K>
class MergeSort {

  /*
   * The superKeyCompare function compares two K arrays in all k dimensions,
   * and uses the sorting or partition coordinate as the most significant dimension.
   *
   * Calling parameters:
   *
   * a - a K *
   * b - a K *
   * p - the most significant dimension
   * dim - the number of dimensions
   *
   * returns a K result of comparing two K arrays
   */
private:
  inline
  static K superKeyCompare(K const* const a,
                           K const* const b,
                           signed_size_t const p,
                           signed_size_t const dim) {
    
    // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
    K diff = a[p] - b[p];
    for (signed_size_t i = 1; diff == 0 && i < dim; ++i) {
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
   * The mergeSortReferenceAscending function recursively subdivides the reference array then
   * merges the elements in ascending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a K ** that represents the array of (x, y, z, w...) coordinates to sort
   * temporary - a K ** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceAscending(K ** const reference,
                                          K ** const temporary,
                                          signed_size_t const low,
                                          signed_size_t const high,
                                          signed_size_t const p,
                                          signed_size_t const dim,
                                          signed_size_t const maximumSubmitDepth,
                                          signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
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
            (superKeyCompare(temporary[i], temporary[j], p, dim) < 0) ? temporary[i++] : temporary[j--];
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

        // Are there sufficient temporary array elements to justify dual-threaded merge?
        if (high - low + 1 > MERGE_CUTOFF)
        {
          // Yes, so compare the results in the temporary array in ascending order with a child thread
          // and merge them into the lower half of the reference array in ascending order.
          auto mergeFuture =
            async(launch::async, [&] {
                                  for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                    reference[k] =
                                      (superKeyCompare(temporary[i], temporary[j], p, dim) <= 0)
                                      ? temporary[i++] : temporary[j--];
                                  }
                                });

          // And simultaneously compare the results in the temporary array in descending order with the
          // current thread and merge them into the upper half of the reference array in ascending order.
          for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
            reference[k] =
              (superKeyCompare(temporary[i], temporary[j], p, dim) > 0) ? temporary[i--] : temporary[j++];
          }

          // Wait for the child thread to finish execution.
          try {
            mergeFuture.get();
          }
          catch (exception const& e) {
            throw runtime_error("\n\ncaught exception for merge future in mergeSortReferenceAscending\n");
          }
        }
        else
        {
          // No, there are insufficient temporary array elements to justify dual-threaded merge,
          // so compare the results in the temporary array in ascending order and merge them into
          // the reference array in ascending order.
          for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
            reference[k] =
              (superKeyCompare(temporary[i], temporary[j], p, dim) < 0) ? temporary[i++] : temporary[j--];
          }
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        K * const tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1], tmp, p, dim) > 0; --j) {
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
   * reference - a K ** that represents the array of (x, y, z, w...) coordinates to sort
   * temporary - a K ** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceDescending(K ** const reference,
                                           K ** const temporary,
                                           signed_size_t const low,
                                           signed_size_t const high,
                                           signed_size_t const p,
                                           signed_size_t const dim,
                                           signed_size_t const maximumSubmitDepth,
                                           signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
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
            (superKeyCompare(temporary[i], temporary[j], p, dim) > 0) ? temporary[i++] : temporary[j--];
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

        // Are there sufficient temporary array elements to justify dual-threaded merge?
        if (high - low + 1 > MERGE_CUTOFF)
        {
          // Yes, so compare the results in the temporary array in ascending order with a child thread
          // and merge them into the lower half of the reference array in descending order.
          auto mergeFuture =
            async(launch::async, [&] {
                                  for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                    reference[k] =
                                      (superKeyCompare(temporary[i], temporary[j], p, dim) >= 0)
                                      ? temporary[i++] : temporary[j--];
                                  }
                                });

          // And simultaneously compare the results in the temporary array in descending order with the
          // current thread and merge them into the upper half of the reference array in descending order.
          for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
            reference[k] =
              (superKeyCompare(temporary[i], temporary[j], p, dim) < 0) ? temporary[i--] : temporary[j++];
          }

          // Wait for the child thread to finish execution.
          try {
            mergeFuture.get();
          }
          catch (exception const& e) {
            throw runtime_error("\n\ncaught exception for merge future in mergeSortReferenceDescending\n");
          }
        }
        else
        {
          // No, there are insufficient temporary array elements to justify dual-threaded merge,
          // so compare the results in the temporary array in ascending order and merge them into
          // the reference array in descending order.
          for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
            reference[k] =
              (superKeyCompare(temporary[i], temporary[j], p, dim) > 0) ? temporary[i++] : temporary[j--];
          }
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
      for (signed_size_t i = low + 1; i <= high; ++i) {
        K * const tmp = reference[i];
        signed_size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1], tmp, p, dim) < 0; --j) {
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
   * reference - a K ** that represents the array of (x, y, z, w...) coordinates to sort
   * temporary - a K ** temporary array into which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryAscending(K ** const reference,
                                          K ** const temporary,
                                          signed_size_t const low,
                                          signed_size_t const high,
                                          signed_size_t const p,
                                          signed_size_t const dim,
                                          signed_size_t const maximumSubmitDepth,
                                          signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
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
            (superKeyCompare(reference[i], reference[j], p, dim) < 0) ? reference[i++] : reference[j--];
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

        // Are there sufficient reference array elements to justify dual-threaded merge?
        if (high - low + 1 > MERGE_CUTOFF)
        {
          // Yes, compare the results in the reference array in ascending order with a child thread
          // and merge them into the lower half of the temporary array in ascending order.
          auto mergeFuture =
            async(launch::async, [&] {
                                  for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                    temporary[k] =
                                      (superKeyCompare(reference[i], reference[j], p, dim) <= 0)
                                      ? reference[i++] : reference[j--];
                                  }
                                });

          // And simultaneously compare the results in the reference array in descending order with the
          // current thread and merge them into the upper half of the temporary array in ascending order.
          for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
            temporary[k] =
              (superKeyCompare(reference[i], reference[j], p, dim) > 0) ? reference[i--] : reference[j++];
          }

          // Wait for the child thread to finish execution.
          try {
            mergeFuture.get();
          }
          catch (exception const& e) {
            throw runtime_error("\n\ncaught exception for merge future in mergeSortTemporaryAscending\n");
          }
        }
        else
        {
          // No, there are insufficient reference array elements to justify dual-threaded merge,
          // so compare the results in the reference array in ascending order and merge them into
          // the temporary array in ascending order.
          for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
            temporary[k] =
              (superKeyCompare(reference[i], reference[j], p, dim) < 0) ? reference[i++] : reference[j--];
          }
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
          if (superKeyCompare(reference[j], temporary[i + 1], p, dim) > 0) {
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
   * reference - a K ** that represents the array of (x, y, z, w...) coordinates to sort
   * temporary - a K ** temporary array into which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryDescending(K ** const reference,
                                           K ** const temporary,
                                           signed_size_t const low,
                                           signed_size_t const high,
                                           signed_size_t const p,
                                           signed_size_t const dim,
                                           signed_size_t const maximumSubmitDepth,
                                           signed_size_t const depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      signed_size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
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
            (superKeyCompare(reference[i], reference[j], p, dim) > 0) ? reference[i++] : reference[j--];
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

        // Are there sufficient reference array elements to justify dual-threaded merge?
        if (high - low + 1 > MERGE_CUTOFF)
        {
          // Yes, so compare the results in the reference array in ascending order with a child thread
          // and merge them into the lower half of the temporary array in descending order.
          auto mergeFuture =
            async(launch::async, [&] {
                                  for (signed_size_t i = low, j = high, k = low; k <= mid; ++k) {
                                    temporary[k] =
                                      (superKeyCompare(reference[i], reference[j], p, dim) >= 0)
                                      ? reference[i++] : reference[j--];
                                  }
                                });

          // And simultaneously compare the results in the reference array in descending order with the
          // current thread and merge them into the upper half of the temporary array in descending order.
          for (signed_size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
            temporary[k] =
              (superKeyCompare(reference[i], reference[j], p, dim) < 0) ? reference[i--] : reference[j++];
          }

          // Wait for the child thread to finish execution.
          try {
            mergeFuture.get();
          }
          catch (exception const& e) {
            throw runtime_error("\n\ncaught exception for merge future in mergeSortTemporaryDescending\n");
          }
        }
        else
        {
          // No, there are insufficient reference array elements to justify dual-threaded merge,
          // sp compare the results in the reference array in ascending order and merge them into
          // the temporary array in descending order.
          for (signed_size_t i = low, j = high, k = low; k <= high; ++k) {
            temporary[k] =
              (superKeyCompare(reference[i], reference[j], p, dim) > 0) ? reference[i++] : reference[j--];
          }
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
          if (superKeyCompare(reference[j], temporary[i + 1], p, dim) < 0) {
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

  friend class KdNode<K>;
  friend class KdTree<K>;
  friend class KdTreeDynamic<K>;
}; // class MergeSort

#endif // KD_TREE_MERGE_SORT_H
