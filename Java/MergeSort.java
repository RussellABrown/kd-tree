/*
 * Copyright (c) 2015, 2019, 2020, 2023, 2025 Russell A. Brown
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

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The {@code MergeSort} class contains multi-threaded merge sort methods.
 * </p>
 */
public class MergeSort
{
    /**
     * <p>
     * The {@code superKeyCompare} method compares two long[] in as few coordinates as possible
     * and uses the sorting or partition coordinate as the most significant coordinate.
     * </p>
     * 
     * @param a - a long[]
     * @param b - a long[]
     * @param p - the most significant dimension
     * @returns a long that represents the result of comparing two super keys
     */
   protected static long superKeyCompare(final long[] a,
                                         final long[] b,
                                         final int p)
    {
        long diff = a[p] - b[p];
        for (int i = 1; diff == 0L && i < a.length; i++) {
            int r = i + p;
            // A fast alternative to the modulus operator for (i + p) < 2 * a.length.
            r = (r < a.length) ? r : r - a.length;
            diff = a[r] - b[r];
        }
        return diff;
    }

    /**
     * <p>
     * The {@code removeDuplicates} method} checks the validity of the merge sort
     * and removes from the reference array all but one of a set of references that
     * reference duplicate tuples.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param p - the index of the most significant coordinate in the super key
     * @returns the address of the last element of the references array following duplicate removal
     */
    protected static int removeDuplicates(final KdNode[] reference,
                                          final int p)
    {
        int end = 0;
        for (int i = 1; i < reference.length; i++) {
            long compare = superKeyCompare(reference[i].tuple, reference[i-1].tuple, p);
            if (compare < 0L) {
                throw new RuntimeException( "merge sort failure: superKeyCompare(ref[" +
                                            Integer.toString(i) + "], ref[" + Integer.toString(i-1) +
                                            "], (" + Integer.toString(p) + ") = " + Long.toString(compare) );
            } else if (compare > 0L) {
                // Keep this element of the reference array by
                // appending it to the kept elements of the array.
                reference[++end] = reference[i];
            } else {
                // Add the values from this element of the reference
                // array to the values set of the duplicate element
                // that is kept, and do not keep this element.
                reference[end].values.addAll(reference[i].values);
            }
        }
        return end;
    }

    /**
     * <p>
     * The following four merge sort methods are adapted from the mergesort function that is shown
     * on p. 166 of Robert Sedgewick's "Algorithms in C++", Addison-Wesley, Reading, MA, 1992.
     * That elegant implementation of the merge sort algorithm eliminates the requirement to test
     * whether the upper and lower halves of an auxiliary array have become exhausted during the
     * merge operation that copies from the auxiliary array to a result array.  This elimination is
     * made possible by inverting the order of the upper half of the auxiliary array and by accessing
     * elements of the upper half of the auxiliary array from highest address to lowest address while
     * accessing elements of the lower half of the auxiliary array from lowest address to highest
     * address.
     * </p>
     * <p>
     * The following four merge sort methods also implement two suggestions from p. 275 of Robert
     * Sedgewick's and Kevin Wayne's "Algorithms 4th Edition", Addison-Wesley, New York, 2011.  The
     * first suggestion is to replace merge sort with insertion sort when the size of the array to
     * sort falls below a threshold.  The second suggestion is to avoid unnecessary copying to the
     * auxiliary array prior to the merge step of the algorithm by implementing two versions of
     * merge sort and by applying some "recursive trickery" to arrange that the required result is
     * returned in the auxiliary array by one version and in the result array by the other version.
     * The following four merge sort methods build upon this suggestion and return their result in
     * either ascending or descending order, as discussed on pp. 173-174 of Robert Sedgewick's
     * Algorithms in C++", Addison-Wesley, Reading, MA, 1992, 
     * </p>
     * <p>
     * During multi-threaded execution, the upper and lower halves of the result array may be filled
     * from the auxiliary array (or vice versa) simultaneously by two threads.  The lower half of the
     * result array is filled by accessing elements of the upper half of the auxiliary array from highest
     * address to lowest address while accessing elements of the lower half of the auxiliary array from
     * lowest address to highest address, as explained above for elimination of the test for exhaustion.
     * The upper half of the result array is filled by addressing elements from the upper half of the
     * auxiliary array from lowest address to highest address while accessing the elements from the lower
     * half of the auxiliary array from highest address to lowest address.  Note: for the upper half
     * of the result array, there is no requirement to test for exhaustion provided that upper half of
     * the result array never comprises more elements than the lower half of the result array.  This
     * provision is satisfied by computing the median address of the result array as shown below for
     * all four merge sort methods. 
     * </p>
     * <p>
     * The {@code mergeSortReferenceAscending} method recursively subdivides the array to be sorted
     * then merges the elements in ascending order and leaves the result in the reference array.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] from which to copy results;
     *                    this array must be as large as the reference array
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    protected static void mergeSortReferenceAscending(final KdNode[] reference,
                                                      final KdNode[] temporary,
                                                      final int low,
                                                      final int high,
                                                      final int p,
                                                      final ExecutorService executor,
                                                      final int maximumSubmitDepth,
                                                      int depth) {

        if (high - low > Constants.INSERTION_SORT_CUTOFF) {

            // Avoid overflow when calculating the median address.
            final int mid = low + ( (high - low) >> 1 );

            // Subdivide the lower half of the tree with a child thread at as many levels of the tree as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.

            // Is a child thread available to subdivide the lower half of the reference array,
            // and are there sufficient KdNode instances to justify spawning a child thread?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth || high - low < Constants.MERGE_CUTOFF) {

                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, low, mid, p, executor, maximumSubmitDepth, depth + 1);

                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Compare the results in the temporary array in ascending order and merge them into
                // the reference array in ascending order.
                for (int i = low, j = high, k = low; k <= high; k++) {
                    reference[k] =
                        (superKeyCompare(temporary[i].tuple, temporary[j].tuple, p) < 0L) ? temporary[i++] : temporary[j--];
                }
                
            } else {

                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the temporary array in ascending order.
                final Future<Void> sortFuture =
                    executor.submit( mergeSortTemporaryAscendingWithThread(reference, temporary,
                                                                            low, mid, p, executor, maximumSubmitDepth, depth + 1) );

                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Then get the result of subdividing the lower half of the reference array with the child thread.
                try {
                    sortFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "sort future exception: " + e.getMessage() );
                }
                
                // Compare the results in the temporary array in ascending order with a child thread
                // and merge them into the lower half of the reference array in ascending order.
                final Future<Void> mergeFuture =
                    executor.submit( mergeResultsAscendingWithThread(reference, temporary, low, high, low, mid, p) );

                // And simultaneously compare the results in the temporary array in descending order with the
                // current thread and merge them into the upper half of the reference array in ascending order.
                for (int i = mid, j = mid + 1, k = high; k > mid; k--) {
                    reference[k] =
                        (superKeyCompare(temporary[i].tuple, temporary[j].tuple, p) > 0L) ? temporary[i--] : temporary[j++];
                }

                // Then get the result of merging into the lower half of the reference array with the child thread.
                try {
                    mergeFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "merge future exception: " + e.getMessage() );
                }
            }

        } else {
            
            // Here is Jon Bentley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
            // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
            for (int i = low + 1; i <= high; i++) {
                KdNode tmp = reference[i];
                int j;
                for (j = i; j > low && superKeyCompare(reference[j-1].tuple, tmp.tuple, p) > 0L; j--) {
                    reference[j] = reference[j-1];
                }
                reference[j] = tmp;
            }
        }
    }
    
    /**
     * <p>
     * The {@code mergeSortReferenceDescending} method recursively subdivides the array to be sorted
     * then merges the elements in descending order and leaves the result in the reference array.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] from which to copy results;
     *                    this array must be as large as the reference array
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static void mergeSortReferenceDescending(final KdNode[] reference,
                                                     final KdNode[] temporary,
                                                     final int low,
                                                     final int high,
                                                     final int p,
                                                     final ExecutorService executor,
                                                     final int maximumSubmitDepth,
                                                     int depth) {

        if (high - low > Constants.INSERTION_SORT_CUTOFF) {

            // Avoid overflow when calculating the median address.
            final int mid = low + ( (high - low) >> 1 );

            // Subdivide the lower half of the tree with a child thread at as many levels of the tree as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.

            // Is a child thread available to subdivide the lower half of the reference array,
            // and are there sufficien KdNode instances to justify spawning a child thread?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth || high - low < Constants.MERGE_CUTOFF) {

                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, low, mid, p, executor, maximumSubmitDepth, depth + 1);

                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Compare the results in the temporary array in ascending order and merge them into
                // the reference array in descending order.
                for (int i = low, j = high, k = low; k <= high; k++) {
                    reference[k] =
                        (superKeyCompare(temporary[i].tuple, temporary[j].tuple, p) > 0L) ? temporary[i++] : temporary[j--];
                }
                
            } else {

                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the temporary array in descending order.
                final Future<Void> sortFuture =
                    executor.submit( mergeSortTemporaryDescendingWithThread(reference, temporary,
                                                                            low, mid, p, executor,
                                                                            maximumSubmitDepth, depth + 1) );

                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Then get the result of subdividing the lower half of the reference array with the child thread.
                try {
                    sortFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "sort future exception: " + e.getMessage() );
                }

                // Compare the results in the temporary array in ascending order with a child thread
                // and merge them into the lower half of the reference array in descending order.
                final Future<Void> mergeFuture =
                    executor.submit( mergeResultsDescendingWithThread(reference, temporary, low, high, low, mid, p) );

                // And simultaneously compare the results in the temporary array in descending order with the
                // current thread and merge them into the upper half of the reference array in descending order.
                for (int i = mid, j = mid + 1, k = high; k > mid; k--) {
                    reference[k] =
                        (superKeyCompare(temporary[i].tuple, temporary[j].tuple, p) < 0L) ? temporary[i--] : temporary[j++];
                }

                // Then get the result of merging into the lower half of the reference array with the child thread.
                try {
                    mergeFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "merge future exception: " + e.getMessage() );
                }
            }

        } else {
            
            // Here is Jon Bentley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
            // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
            for (int i = low + 1; i <= high; i++) {
                KdNode tmp = reference[i];
                int j;
                for (j = i; j > low && superKeyCompare(reference[j-1].tuple, tmp.tuple, p) < 0; j--) {
                    reference[j] = reference[j-1];
                }
                reference[j] = tmp;
            }
        }
    }

    /**
     * <p>
     * The {@code mergeSortTemporaryAscending} method recursively subdivides the array to be sorted
     * then merges the elements in ascending order and leaves the result in the temporary array.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] into which to copy results;
     *                    this array must be as large as the reference array
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static void mergeSortTemporaryAscending(final KdNode[] reference,
                                                    final KdNode[] temporary,
                                                    final int low,
                                                    final int high,
                                                    final int p,
                                                    final ExecutorService executor,
                                                    final int maximumSubmitDepth,
                                                    int depth) {

        if (high - low > Constants.INSERTION_SORT_CUTOFF) {

            // Avoid overflow when calculating the median address.
            final int mid = low + ( (high - low) >> 1 );

            // Subdivide the lower half of the tree with a child thread at as many levels of the tree as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.

            // Is a child thread available to subdivide the lower half of the reference array,
            // and are there sufficien KdNode instances to justify spawning a child thread?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth || high - low < Constants.MERGE_CUTOFF) {

                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, low, mid, p, executor, maximumSubmitDepth, depth + 1);

                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Compare the results in the reference array in ascending order and merge them into
                // the temporary array in ascending order.
                for (int i = low, j = high, k = low; k <= high; k++) {
                    temporary[k] =
                        (superKeyCompare(reference[i].tuple, reference[j].tuple, p) < 0L) ? reference[i++] : reference[j--];
                }

            
            } else {

                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the reference array in ascending order.
                final Future<Void> sortFuture =
                    executor.submit( mergeSortReferenceAscendingWithThread(reference, temporary,
                                                                            low, mid, p, executor,
                                                                            maximumSubmitDepth, depth + 1) );

                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Then get the result of subdividing the lower half of the reference array with the child thread.
                try {
                    sortFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "sort future exception: " + e.getMessage() );
                }

                // Compare the results in the reference array in ascending order with a child thread
                // and merge them into the lower half of the temporary array in ascending order.
                final Future<Void> mergeFuture =
                    executor.submit( mergeResultsAscendingWithThread(temporary, reference, low, high, low, mid, p) );

                // And simultaneously compare the results in the reference array in descending order with the
                // current thread and merge them into the upper half of the temporary array in ascending order.
                for (int i = mid, j = mid + 1, k = high; k > mid; k--) {
                    temporary[k] =
                        (superKeyCompare(reference[i].tuple, reference[j].tuple, p) > 0L) ? reference[i--] : reference[j++];
                }

                // Then get the result of merging into the lower half of the temporary array with the child thread.
                try {
                    mergeFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "merge future exception: " + e.getMessage() );
                }
            }

        } else {
            
            // This implementation of insertion sort leaves the result in the temporary array in ascending order.
            temporary[high] = reference[high];
            int i, j;
            for (j = high - 1; j >= low; j--) {
                for (i = j; i < high; i++) {
                    if(superKeyCompare(reference[j].tuple, temporary[i + 1].tuple, p) > 0L) {
                        temporary[i] = temporary[i + 1];
                    } else {
                        break;
                    }
                }
                temporary[i] = reference[j];
            }
        }
    }
    
    /**
     * <p>
     * The {@code mergeSortTemporaryDescending} method recursively subdivides the array to be sorted
     * then merges the elements in descending order and leaves the result in the temporary array.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch array into which to copy results;
     *                    this array must be as large as the reference array
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static void mergeSortTemporaryDescending(final KdNode[] reference,
                                                     final KdNode[] temporary,
                                                     final int low,
                                                     final int high,
                                                     final int p,
                                                     final ExecutorService executor,
                                                     final int maximumSubmitDepth,
                                                     int depth) {

        if (high - low > Constants.INSERTION_SORT_CUTOFF) {

            // Avoid overflow when calculating the median address.
            final int mid = low + ( (high - low) >> 1 );

            // Subdivide the lower half of the tree with a child thread at as many levels of the tree as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.

            // Is a child thread available to subdivide the lower half of the reference array,
            // and are there sufficien KdNode instances to justify spawning a child thread?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth || high - low < Constants.MERGE_CUTOFF) {

                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, low, mid, p, executor, maximumSubmitDepth, depth + 1);

                // Then recursively subdivide the upper half of the reference array with the current thread
                // thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Compare the results in the reference array in ascending order and merge them into
                // the temporary array in descending order.
                for (int i = low, j = high, k = low; k <= high; k++) {
                    temporary[k] =
                        (superKeyCompare(reference[i].tuple, reference[j].tuple, p) > 0L) ? reference[i++] : reference[j--];
                }

                
            } else {

                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the reference array in descending order.
                final Future<Void> sortFuture =
                    executor.submit( mergeSortReferenceDescendingWithThread(reference, temporary,
                                                                            low, mid, p, executor,
                                                                            maximumSubmitDepth, depth + 1) );

                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, executor, maximumSubmitDepth, depth + 1);

                // Then get the result of subdividing the lower half of the reference array with the child thread.
                try {
                    sortFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "sort future exception: " + e.getMessage() );
                }

                // Compare the results in the reference array in ascending order with a child thread
                // and merge them into the lower half of the temporary array in descending order.
                final Future<Void> mergeFuture =
                    executor.submit( mergeResultsDescendingWithThread(temporary, reference, low, high, low, mid, p) );

                // And simultaneously compare the results in the reference array in descending order with the
                // current thread and merge them into the upper half of the temporary array in descending order.
                for (int i = mid, j = mid + 1, k = high; k > mid; k--) {
                    temporary[k] =
                        (superKeyCompare(reference[i].tuple, reference[j].tuple, p) < 0L) ? reference[i--] : reference[j++];
                }

                // Then get the result of merging into the lower half of the temporary array with the child thread.
                try {
                    mergeFuture.get();
                } catch (Exception e) {
                    throw new RuntimeException( "merge future exception: " + e.getMessage() );
                }
            }

        } else {
            
            // This implementation of insertion sort leaves the result in the temporary array in descending order.
            temporary[high] = reference[high];
            int i, j;
            for (j = high - 1; j >= low; j--) {
                for (i = j; i < high; i++) {
                    if (superKeyCompare(reference[j].tuple, temporary[i + 1].tuple, p) < 0L) {
                        temporary[i] = temporary[i + 1];
                    } else {
                        break;
                    }
                }
                temporary[i] = reference[j];
            }
        }
    }
    
    /**
     * <p> The {@code mergeResultsAscending} method compares the results in the source array in order of
     * ascending address and merges them into the destination array in order of ascending value.
     * </p>
     * 
     * @param destination - a KdNode[] from which to merge results
     * @param source - a KdNode[] into which to merge results
     * @param iStart - the initial value of the i-index
     * @param jStart - the initial value of the j-index
     * @param kStart - the initial value of the k-index
     * @param kEnd - the final value of the k-index
     * @param p - the sorting partition (x, y, z, w...)
     */
    private static void mergeResultsAscending(final KdNode[] destination,
                                              final KdNode[] source,
                                              final int iStart,
                                              final int jStart,
                                              final int kStart,
                                              final int kEnd,
                                              final int p) {
        
        for (int i = iStart, j = jStart, k = kStart; k <= kEnd; k++) {
            destination[k] = (superKeyCompare(source[i].tuple, source[j].tuple, p) <= 0) ? source[i++] : source[j--];
        }
    }

    /**
     * <p> The {@code mergeResultsDescending} method compares the results in the source array in order of
     * ascending address and merges them into the destination array in order of descending value.
     * </p>
     * 
     * @param destination - a KdNode[] from which to merge results
     * @param source - a KdNode[] into which to merge results
     * @param iStart - the initial value of the i-index
     * @param jStart - the initial value of the j-index
     * @param kStart - the initial value of the k-index
     * @param kEnd - the final value of the k-index
     * @param p - the sorting partition (x, y, z, w...)
     */
    private static void mergeResultsDescending(final KdNode[] destination,
                                               final KdNode[] source,
                                               final int iStart,
                                               final int jStart,
                                               final int kStart,
                                               final int kEnd,
                                               final int p) {
        
        for (int i = iStart, j = jStart, k = kStart; k <= kEnd; k++) {
            destination[k] = (superKeyCompare(source[i].tuple, source[j].tuple, p) >= 0) ? source[i++] : source[j--];
        }
    }

    /**
     * <p>
     * The {@code mergeResultsAscendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeResultsAscending mergeResultsAscending} method.
     * </p>
     * 
     * @param destination - a KdNode[] from which to merge results
     * @param source - a KdNode[] into which to merge results
     * @param iStart - the initial value of the i-index
     * @param jStart - the initial value of the j-index
     * @param kStart - the initial value of the k-index
     * @param kEnd - the final value of the k-index
     * @param p - the sorting partition (x, y, z, w...)
     */
    private static Callable<Void> mergeResultsAscendingWithThread(final KdNode[] destination,
                                                                  final KdNode[] source,
                                                                  final int iStart,
                                                                  final int jStart,
                                                                  final int kStart,
                                                                  final int kEnd,
                                                                  final int p) {
        
        return new Callable<Void>() {
            @Override
            public Void call() {
                mergeResultsAscending(destination, source, iStart, jStart, kStart, kEnd, p);
                return null;
            }
        };
    }

    /**
     * <p>
     * The {@code mergeResultsDescendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeResultsDescending mergeResultsDescending} method.
     * </p>
     * 
     * @param destination - a KdNode[] from which to merge results
     * @param source - a KdNode[] into which to merge results
     * @param iStart - the initial value of the i-index
     * @param jStart - the initial value of the j-index
     * @param kStart - the initial value of the k-index
     * @param kEnd - the final value of the k-index
     * @param p - the sorting partition (x, y, z, w...)
     */
    private static Callable<Void> mergeResultsDescendingWithThread(final KdNode[] destination,
                                                                   final KdNode[] source,
                                                                   final int iStart,
                                                                   final int jStart,
                                                                   final int kStart,
                                                                   final int kEnd,
                                                                   final int p) {
        
        return new Callable<Void>() {
            @Override
            public Void call() {
                mergeResultsDescending(destination, source, iStart, jStart, kStart, kEnd, p);
                return null;
            }
        };
    }

    /**
     * <p>
     * The {@code mergeSortReferenceAscendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeSortReferenceAscending mergeSortReferenceAscending} method.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] from which to copy results;
     *                    this array must be as large as the reference array.
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static Callable<Void> mergeSortReferenceAscendingWithThread(final KdNode[] reference,
                                                                        final KdNode[] temporary,
                                                                        final int low,
                                                                        final int high,
                                                                        final int p,
                                                                        final ExecutorService executor,
                                                                        final int maximumSubmitDepth,
                                                                        final int depth) {
        
        return new Callable<Void>() {
            @Override
            public Void call() {
                mergeSortReferenceAscending(reference, temporary, low, high, p, executor, maximumSubmitDepth, depth);
                return null;
            }
            };
    }

    /**
     * <p>
     * The {@code mergeSortReferenceDescendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeSortReferenceDescending mergeSortReferenceDescending} method.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] from which to copy results;
     *                    this array must be as large as the reference array.
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static Callable<Void> mergeSortReferenceDescendingWithThread(final KdNode[] reference,
                                                                         final KdNode[] temporary,
                                                                         final int low,
                                                                         final int high,
                                                                         final int p,
                                                                         final ExecutorService executor,
                                                                         final int maximumSubmitDepth,
                                                                         final int depth) {
        
        return new Callable<Void>() {
            @Override
            public Void call() {
                mergeSortReferenceDescending(reference, temporary, low, high, p, executor, maximumSubmitDepth, depth);
                return null;
            }
        };
    }

    /**
     * <p>
     * The {@code mergeSortTemporaryAscendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeSortTemporaryAscending mergeSortTemporaryAscending} method.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] into which to copy results;
     *                    this array must be as large as the reference array.
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference array
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static Callable<Void> mergeSortTemporaryAscendingWithThread(final KdNode[] reference,
                                                                        final KdNode[] temporary,
                                                                        final int low,
                                                                        final int high,
                                                                        final int p,
                                                                        final ExecutorService executor,
                                                                        final int maximumSubmitDepth,
                                                                        final int depth) {
        
        return new Callable<Void>() {
            @Override
            public Void call() {
                mergeSortTemporaryAscending(reference, temporary, low, high, p, executor, maximumSubmitDepth, depth);
                return null;
            }
        };
    }

    /**
     * <p>
     * The {@code mergeSortTemporaryDescendingWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#mergeSortTemporaryDescending mergeSortTemporaryDescending} method.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] into which to copy results;
     *                    this array must be as large as the reference array.
     * @param low - the start index of the region of the reference array
     * @param high - the high index of the region of the reference
     * @param p - the sorting partition (x, y, z, w...)
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param depth - the depth of subdivision
     */
    private static Callable<Void> mergeSortTemporaryDescendingWithThread(final KdNode[] reference,
                                                                         final KdNode[] temporary,
                                                                         final int low,
                                                                         final int high,
                                                                         final int p,
                                                                         final ExecutorService executor,
                                                                         final int maximumSubmitDepth,
                                                                         final int depth) {
        
        return new Callable<Void>() {
            @Override
                public Void call() {
                    mergeSortTemporaryDescending(reference, temporary, low, high, p, executor, maximumSubmitDepth, depth);
                    return null;
                }
        };
    }
} // class MergeSort
