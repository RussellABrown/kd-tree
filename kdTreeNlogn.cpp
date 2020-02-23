/*
 * Copyright (c) 2015, Russell A. Brown
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
 * Gnu g++ compilation options are: -lm -O3 (compile using g++ not gcc)
 *
 * Optional compilation options are:
 *
 * -DMACH - Use a Mach equivalent to the clock_gettime(CLOCK_REALTIME, &time) function.
 * -DDUAL_THREAD_MEDIAN - Calculate the medians with two threads.
 * -DDUAL_THREAD_INDEX - Find the index of the median of medians with two threads.
 * -DBIDIRECTIONAL_PARTITION - Partition an array about the median of medians proceeding
 *                             from both ends of the array instead of only the beginning.
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

using namespace std;

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
struct timespec getTime(void) {
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return time;
}
#endif

/* One node of a k-d tree */
template <typename T>
class KdNode
{
private:
    T const* tuple;
    const KdNode *ltChild, *gtChild;

public:
    KdNode() {
        this->ltChild = this->gtChild = NULL;
    }

public:
    T const* getTuple() const {
        return this->tuple;
    }
    
    /*
     * The getKdNode function gets a KdNode from the kdNodes array and assigns the tuple field.
     *
     * calling parameters:
     *
     * reference - a T** that represents one of the reference arrays
     * kdNodes - a vector<KdNode*> that contains pre-allocated KdNodes
     * k - the index into both the reference and the kdNodes arrays
     *
     * returns: a KdNode* to the KdNode to which the tuple has been assigned
     */
private:
    inline
    static KdNode* getKdNode(T** reference, vector<KdNode*>& kdNodes, long const k) {
        KdNode* kdNode = kdNodes[k];
        kdNode->tuple = reference[k];
        return kdNode;
    }
    
    /*
     * The initializeReferences function initializes a reference array by creating
     * references into the coordinates array.
     *
     * calling parameters:
     *
     * coordinates - a vector<T*> of pointers to each of the (x, y, z, w...) tuples
     * reference - a T** that represents one of the reference arrays
     */
private:
    inline
    static void initializeReference(vector<T*>& coordinates, T** reference) {
        for (size_t j = 0; j < coordinates.size(); ++j) {
            reference[j] = coordinates[j];
        }
    }
    
    /*
     * The superKeyCompare function compares two T arrays in all k dimensions,
     * and uses the sorting or partition coordinate as the most significant dimension.
     *
     * calling parameters:
     *
     * a - a T*
     * b - a T*
     * p - the most significant dimension
     * dim - the number of dimensions
     *
     * returns: a long result of comparing two T arrays
     */
private:
    inline
    static T superKeyCompare(T const* a, T const*b, long const p, long const dim) {
        // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
        T diff = a[p] - b[p];
        for (long i = 1; diff == 0 && i < dim; i++) {
            long r = i + p;
            // A fast alternative to the modulus operator for (i + p) < 2 * dim.
            r = (r < dim) ? r : r - dim;
            diff = a[r] - b[r];
        }
        return diff;
    }

    /* A cutoff for switching from merge sort to insertion sort */
    const static long INSERTION_SORT_CUTOFF = 15;

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
     * calling parameters:
     *
     * reference - a T** that represents the array of (x, y, z, w...) coordinates to sort
     * temporary - a T** temporary array from which to copy sorted results;
     *             this array must be as large as the reference array
     * low - the start index of the region of the reference array to sort
     * high - the end index of the region of the reference array to sort
     * p - the sorting partition (x, y, z, w...)
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the tree depth
     */
private:
    static void mergeSortReferenceAscending(T** reference, T** temporary,
                                            long const low, long const high, long const p, long const dim,
                                            long const maximumSubmitDepth, long const depth) {
                                                     
        if (high - low > INSERTION_SORT_CUTOFF) {
                 
            // Avoid overflow when calculating the median.
            long const mid = low + ( (high - low) >> 1 );
                 
            // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
            // Is a child thread available to subdivide the lower half of the reference array?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                     
                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Compare the results in the temporary array in ascending order and merge them into
                // the reference array in ascending order.
                for (long i = low, j = high, k = low; k <= high; ++k) {
                    reference[k] =
                        (superKeyCompare(temporary[i], temporary[j], p, dim) < 0) ? temporary[i++] : temporary[j--];
                }
                     
            } else {
                     
                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the temporary array in ascending order.
                future<void> sortFuture = async(launch::async, mergeSortTemporaryAscending, reference, temporary,
                                                low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Wait for the child thread to finish execution.
                try {
                    sortFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
                     
                // Compare the results in the temporary array in ascending order with a child thread
                // and merge them into the lower half of the reference array in ascending order.
                future<void> mergeFuture = async(launch::async, [&] {
                    for (long i = low, j = high, k = low; k <= mid; ++k) {
                        reference[k] =
                            (superKeyCompare(temporary[i], temporary[j], p, dim) <= 0) ? temporary[i++] : temporary[j--];
                    }
                });
                     
                // And simultaneously compare the results in the temporary array in descending order with the
                // current thread and merge them into the upper half of the reference array in ascending order.
                for (long i = mid, j = mid + 1, k = high; k > mid; --k) {
                    reference[k] =
                        (superKeyCompare(temporary[i], temporary[j], p, dim) > 0) ? temporary[i--] : temporary[j++];
                }
                     
                // Wait for the child thread to finish execution.
                try {
                    mergeFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }
                 
        } else {
                 
            // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
            // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
            for (long i = low + 1; i <= high; ++i) {
                T* tmp = reference[i];
                long j;
                for (j = i; j > low && superKeyCompare(reference[j-1], tmp, p, dim) > 0; j--) {
                    reference[j] = reference[j-1];
                }
                reference[j] = tmp;
            }
        }
    }

    /*
     * The mergeSortReferenceDecending function recursively subdivides the reference array then
     * merges the elements in descending order and leaves the result in the reference array.
     *
     * calling parameters:
     *
     * reference - a T** that represents the array of (x, y, z, w...) coordinates to sort
     * temporary - a T** temporary array from which to copy sorted results;
     *             this array must be as large as the reference array
     * low - the start index of the region of the reference array to sort
     * high - the end index of the region of the reference array to sort
     * p - the sorting partition (x, y, z, w...)
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the tree depth
     */
private:
    static void mergeSortReferenceDescending(T** reference, T** temporary,
                                             long const low, long const high, long const p, long const dim,
                                             long const maximumSubmitDepth, long const depth) {
             
        if (high - low > INSERTION_SORT_CUTOFF) {
                 
            // Avoid overflow when calculating the median.
            long const mid = low + ( (high - low) >> 1 );
                 
            // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
            // Is a child thread available to subdivide the lower half of the reference array?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                     
                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the temporary array in descending order.
                mergeSortTemporaryDescending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Compare the results in the temporary array in ascending order and merge them into
                // the reference array in descending order.
                for (long i = low, j = high, k = low; k <= high; ++k) {
                    reference[k] =
                        (superKeyCompare(temporary[i], temporary[j], p, dim) > 0) ? temporary[i++] : temporary[j--];
                }
                     
            } else {
                     
                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the temporary array in descending order.
                future<void> sortFuture = async(launch::async, mergeSortTemporaryDescending, reference, temporary,
                                                low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the temporary array in ascending order.
                mergeSortTemporaryAscending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Wait for the child thread to finish execution.
                try {
                    sortFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }

                // Compare the results in the temporary array in ascending order with a child thread
                // and merge them into the lower half of the reference array in descending order.
                future<void> mergeFuture = async(launch::async, [&] {
                    for (long i = low, j = high, k = low; k <= mid; ++k) {
                        reference[k] =
                            (superKeyCompare(temporary[i], temporary[j], p, dim) >= 0) ? temporary[i++] : temporary[j--];
                    }
                });
                     
                // And simultaneously compare the results in the temporary array in descending order with the
                // current thread and merge them into the upper half of the reference array in descending order.
                for (long i = mid, j = mid + 1, k = high; k > mid; --k) {
                    reference[k] =
                        (superKeyCompare(temporary[i], temporary[j], p, dim) < 0) ? temporary[i--] : temporary[j++];
                }
                     
                // Wait for the child thread to finish execution.
                try {
                    mergeFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }
                 
        } else {
                 
            // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
            // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
            for (long i = low + 1; i <= high; ++i) {
                T* tmp = reference[i];
                long j;
                for (j = i; j > low && superKeyCompare(reference[j-1], tmp, p, dim) < 0; j--) {
                    reference[j] = reference[j-1];
                }
                reference[j] = tmp;
            }
        }
    }
                                                 
    /*
     * The mergeSortTemporaryAscending function recursively subdivides the reference array then
     * merges the elements in ascending order and leaves the result in the temporary array.
     *
     * calling parameters:
     *
     * reference - a T** that represents the array of (x, y, z, w...) coordinates to sort
     * temporary - a T** temporary array into which to copy sorted results;
     *             this array must be as large as the reference array
     * low - the start index of the region of the reference array to sort
     * high - the end index of the region of the reference array to sort
     * p - the sorting partition (x, y, z, w...)
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the tree depth
     */
private:
    static void mergeSortTemporaryAscending(T** reference, T** temporary,
                                            long const low, long const high, long const p, long const dim,
                                            long const maximumSubmitDepth, long const depth) {
                                                     
        if (high - low > INSERTION_SORT_CUTOFF) {
                 
            // Avoid overflow when calculating the median.
            long const mid = low + ( (high - low) >> 1 );
                 
            // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
            // Is a child thread available to subdivide the lower half of the reference array?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                     
                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Compare the results in the reference array in ascending order and merge them into
                // the temporary array in ascending order.
                for (long i = low, j = high, k = low; k <= high; ++k) {
                    temporary[k] =
                        (superKeyCompare(reference[i], reference[j], p, dim) < 0) ? reference[i++] : reference[j--];
                }
                     
            } else {
                     
                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the reference array in ascending order.
                future<void> sortFuture = async(launch::async, mergeSortReferenceAscending, reference, temporary,
                                                low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Wait for the child thread to finish execution.
                try {
                    sortFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
                     
                // Compare the results in the reference array in ascending order with a child thread
                // and merge them into the lower half of the temporary array in ascending order.
                future<void> mergeFuture = async(launch::async, [&] {
                    for (long i = low, j = high, k = low; k <= mid; ++k) {
                        temporary[k] =
                            (superKeyCompare(reference[i], reference[j], p, dim) <= 0) ? reference[i++] : reference[j--];
                    }
                });

                // And simultaneously compare the results in the reference array in descending order with the
                // current thread and merge them into the upper half of the temporary array in ascending order.
                for (long i = mid, j = mid + 1, k = high; k > mid; --k) {
                    temporary[k] =
                        (superKeyCompare(reference[i], reference[j], p, dim) > 0) ? reference[i--] : reference[j++];
                }
                     
                // Wait for the child thread to finish execution.
                try {
                    mergeFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }
                 
        } else {
                 
            // Here is John Robinson's implementation of insertion sort that sorts in ascending order
            // and leaves the result in the temporary array.
            temporary[high] = reference[high];
            long i, j;
            for (j = high - 1; j >= low; j--){
                for(i = j; i < high; ++i) {
                    if(superKeyCompare(reference[j], temporary[i + 1], p, dim) > 0) {
                        temporary[i] = temporary[i + 1];
                    } else {
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
     * calling parameters:
     *
     * reference - a T** that represents the array of (x, y, z, w...) coordinates to sort
     * temporary - a T** temporary array into which to copy sorted results;
     *             this array must be as large as the reference array
     * low - the start index of the region of the reference array to sort
     * high - the end index of the region of the reference array to sort
     * p - the sorting partition (x, y, z, w...)
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the tree depth
     */
private:
    static void mergeSortTemporaryDescending(T** reference, T** temporary,
                                             long const low, long const high, long const p, long const dim,
                                             long const maximumSubmitDepth, long const depth) {
             
        if (high - low > INSERTION_SORT_CUTOFF) {
                 
            // Avoid overflow when calculating the median.
            long const mid = low + ( (high - low) >> 1 );
                 
            // Subdivide the lower half of the array with a child thread at as many levels of subdivision as possible.
            // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
            // Is a child thread available to subdivide the lower half of the reference array?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                     
                // No, recursively subdivide the lower half of the reference array with the current
                // thread and return the result in the reference array in descending order.
                mergeSortReferenceDescending(reference, temporary, low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // Then recursively subdivide the upper half of the reference array with the current
                // thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Compare the results in the reference array in ascending order and merge them into
                // the temporary array in descending order.
                for (long i = low, j = high, k = low; k <= high; ++k) {
                    temporary[k] =
                        (superKeyCompare(reference[i], reference[j], p, dim) > 0) ? reference[i++] : reference[j--];
                }
                     
            } else {
                     
                // Yes, a child thread is available, so recursively subdivide the lower half of the reference
                // array with a child thread and return the result in the reference array in descending order.
                future<void> sortFuture = async(launch::async, mergeSortReferenceDescending, reference, temporary,
                                                low, mid, p, dim, maximumSubmitDepth, depth+1);
                     
                // And simultaneously, recursively subdivide the upper half of the reference array with
                // the current thread and return the result in the reference array in ascending order.
                mergeSortReferenceAscending(reference, temporary, mid+1, high, p, dim, maximumSubmitDepth, depth+1);
                     
                // Wait for the child thread to finish execution.
                try {
                    sortFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
                     
                // Compare the results in the reference array in ascending order with a child thread
                // and merge them into the lower half of the temporary array in descending order.
                future<void> mergeFuture = async(launch::async, [&] {
                    for (long i = low, j = high, k = low; k <= mid; ++k) {
                        temporary[k] =
                            (superKeyCompare(reference[i], reference[j], p, dim) >= 0) ? reference[i++] : reference[j--];
                    }
                });

                // And simultaneously compare the results in the reference array in descending order with the
                // current thread and merge them into the upper half of the temporary array in descending order.
                for (long i = mid, j = mid + 1, k = high; k > mid; --k) {
                    temporary[k] =
                        (superKeyCompare(reference[i], reference[j], p, dim) < 0) ? reference[i--] : reference[j++];
                }
                     
                // Wait for the child thread to finish execution.
                try {
                    mergeFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }
                 
        } else {
                 
            // Here is John Robinson's implementation of insertion sort that sorts in descending order
            // and leaves the result in the temporary array.
            temporary[high] = reference[high];
            long i, j;
            for (j = high - 1; j >= low; j--){
                for(i = j; i < high; ++i) {
                    if(superKeyCompare(reference[j], temporary[i + 1], p, dim) < 0) {
                        temporary[i] = temporary[i + 1];
                    } else {
                        break;
                    }
                }
                temporary[i] = reference[j];
            }
        }
    }

    /*
     * The removeDuplicates function checks the validity of the merge sort and
     * removes duplicates from a reference array.
     *
     * calling parameters:
     *
     * reference - a T** that represents one of the reference arrays
     * i - the leading dimension for the super key
     * dim - the number of dimensions
     *
     * returns: the end index of the reference array following removal of duplicate elements
     */
private:
    inline
    static long removeDuplicates(T** reference, long const i, long const dim, long const size) {
        long end = 0;
        for (long j = 1; j < size; ++j) {
            long compare = superKeyCompare(reference[j], reference[j-1], i, dim);
            if (compare < 0) {
                cout << "merge sort failure: superKeyCompare(ref[" << j << "], ref["
                << j-1 << "], (" << i << ") = " << compare  << endl;
                exit(1);
            } else if (compare > 0) {
                reference[++end] = reference[j];
            }
        }
        return end;
    }
    
    /*
     * The swap function swaps two array elements.
     *
     * calling parameters:
     *
     * a - array of references to the (x,y,z,w...) coordinates
     * i - the index of the first element
     * j - the index of the second element
     */
private:
    inline
    static void swap(T** a, const long i, const long j) {
        T* t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
		
    /*
     * The following select_j_k functions select the jth of k items.  Adapted
     * from Chapter 4, "Linear Orderings", of Alexander Stepanov's and
     * Paul McJones' "Elements of Programming", Addison-Wesley, New York, 2009.
     *
     * calling parameters:
     *
     * a through e - arrays to include in the selection
     * p - the sorting partition (x, y, z, w...)
     * dim - the number of dimensions
     *
     * returns - a T* that represents the selected array
     */
private:
    inline
    static T* select_0_2(T* a,
                         T* b,
                         const long p,
                         const long dim) {
        if (superKeyCompare(a, b, p, dim) < 0) {
            // a < b
            return a;
        } else {
            // b < a
            return b;
        }
    }

private:
    inline
    static T* select_1_2(T* a,
                         T* b,
                         const long p,
                         const long dim) {
        if (superKeyCompare(a, b, p, dim) < 0) {
            // a < b
            return b;
        } else {
            // b < a
            return a;
        }
    }
		
private:
    inline
    static T* select_1_3_ab(T* a,
                            T* b,
                            T* c,
                            const long p,
                            const long dim) {
        if (superKeyCompare(b, c, p, dim) < 0) {
            // a < b < c
            return b;
        } else {
            // a ? c < b
            return select_1_2(a, c, p, dim);
        }
    }

private:
    inline
    static T* select_1_3(T* a,
                         T* b,
                         T* c,
                         const long p,
                         const long dim) {
        if (superKeyCompare(a, b, p, dim) < 0) {
            // a < b
            return select_1_3_ab(a, b, c, p, dim);
        } else {
            // b < a
            return select_1_3_ab(b, a, c, p, dim);
        }
    }

private:
    inline
    static T* select_1_4_ab_cd(T* a,
                               T* b,
                               T* c,
                               T* d,
                               const long p,
                               const long dim) {
        if (superKeyCompare(c, a, p, dim) < 0) {
            // c < a < b && a ? d so c is eliminated and a ? d
            return select_0_2(a, d, p, dim);
        } else {
            // a < b ? c < d so a is eliminated and b ? c
            return select_0_2(b, c, p, dim);
        }
    }
		
private:
    inline
    static T* select_1_4_ab(T* a,
                            T* b,
                            T* c,
                            T* d,
                            const long p,
                            const long dim) {
        if (superKeyCompare(c, d, p, dim) < 0) {
            // a < b && c < d
            return select_1_4_ab_cd(a, b, c, d, p, dim);
        } else {
            // a < b && d < c
            return select_1_4_ab_cd(a, b, d, c, p, dim);
        }
    }
		
private:
    inline
    static T* select_1_4(T* a,
                         T* b,
                         T* c,
                         T* d,
                            const long p,
                            const long dim) {
        if (superKeyCompare(a, b, p, dim) < 0) {
            // a < b
            return select_1_4_ab(a, b, c, d, p, dim);
        } else {
            // b < a
            return select_1_4_ab(b, a, c, d, p, dim);
        }
    }
		
private:
    inline
    static T* select_2_5_ab_cd(T* a,
                               T* b,
                               T* c,
                               T* d,
                               T* e,
                               const long p,
                               const long dim) {
        if (superKeyCompare(c, a, p, dim) < 0) {
            // c < a < b && c < d ? e so c is eliminated and a < b && d ? e
            return select_1_4_ab(a, b, d, e, p, dim);
        } else {
            // a < b ? c && c < d ? e && b ? e so a is eliminated and c < d && b ? e
            return select_1_4_ab(c, d, b, e, p, dim);
        }
    }
		
private:
    inline
    static T* select_2_5_ab(T* a,
                            T* b,
                            T* c,
                            T* d,
                            T* e,
                            const long p,
                            const long dim) {
        if (superKeyCompare(c, d, p, dim) < 0) {
            // a < b && c < d
            return select_2_5_ab_cd(a, b, c, d, e, p, dim);
        } else {
            // a < b && d < c
            return select_2_5_ab_cd(a, b, d, c, e, p, dim);
        }
    }
		
private:
    inline
    static T* select_2_5(T* a,
                         T* b,
                         T* c,
                         T* d,
                         T* e,
                         const long p,
                         const long dim) {
        if (superKeyCompare(a, b, p, dim) < 0) {
            // a < b
            return select_2_5_ab(a, b, c, d, e, p, dim);
        } else {
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
    static long partition(T** a, const long start, const long n, const long size, const long k,
                          T** medians, const long first, const long p, const long dim,
                          const bool twoThreads) {
			
        if (n <=0 || n > size) {
            cout << "Error in n = " << n << "  size = " << size << endl;
        }
        if (k <= 0 || k > n) {
            cout << "Error in k = " << k << endl;
        }
        if (start + n > size) {
            cout << "Error in start = " << start << "  n = " << n << "  size = " << size << endl;
        }
			
        // This trivial case terminates recursion.
        if ( n == 1 && k == 1 ) {
            return start;
        }

        // Use insertion sort instead of the median of medians algorithm for a small number of elements,
        // via Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
        // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the array a.
        const long MEDIAN_OF_MEDIANS_CUTOFF = 15;
        if (n <= MEDIAN_OF_MEDIANS_CUTOFF) {
            for (long i = start + 1; i <= start + n - 1; ++i) {
                T* tmp = a[i];
                long j;
                for (j = i; j > start && superKeyCompare(a[j-1], tmp, p, dim) > 0; j--) {
                    a[j] = a[j-1];
                }
                a[j] = tmp;
            }
            return start + k - 1;
        }
		    
        // Otherwise, determine how many medians to find.  Round down to count
        // only groups that comprise fully GROUP_SIZE elements.  Any remaining
        // group of elements that doesn't comprise GROUP_SIZE elements will
        // be processed after the following 'for' loop.
        const long GROUP_SIZE = 5; // Must be 5 due to select_2_5 function below.
        long m = n / GROUP_SIZE;
        long startOfGroup;
 
#ifdef DUAL_THREAD_MEDIAN
        // Is more than one thread available to calculate the medians and are
        // there sufficient medians to justify multi-threaded processing?
        //
        // NOTE, however, that NO value of MEDIAN_CUTOFF appears to improve
        // the performance of two threads relative to that of one thread.
        // Hence, the cost of spawning a child thread appears to exceed any
        // improvement in the performance of calculating the medians that
        // may be achieved via two threads.
        const long MEDIAN_CUTOFF = 16384;
        if (twoThreads && m > MEDIAN_CUTOFF) {

            // Yes, calculate the relative index of the middle median.
            long mid = (m + 1) >> 1;
            startOfGroup = mid * GROUP_SIZE;

            // Calculate the lower set of medians with a child thread.
            future<void> medianFuture = async(launch::async, [&] {
                for (long firstOfGroup = 0, i = 0; i < mid; ++i) {
                        
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
            for (long i = mid; i < m; ++i) {
                        
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
            } catch (const exception& e) {
                cout << "caught exception " << e.what() << endl;
            }
        } else
#endif
        {
            // No, only one thread is available, so calculate all medians with the current thread.
            startOfGroup = 0;
            for (long i = 0; i < m; ++i) {
                        
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
        long remainingElements = n - startOfGroup;
        if ( remainingElements < 0 || remainingElements >= GROUP_SIZE ) {
            cout << "Error: incorrect group calculation";
        }
	    	
        // Find the median of any remaining elements via select_j_k.
        switch (remainingElements) {
            case 0:
                break;
            case 1:
                medians[first + m] = a[start + startOfGroup];
                m++;
                break;
            case 2:
                medians[first + m] = select_0_2(a[start + startOfGroup],
                                                a[start + startOfGroup + 1],
                                                p,
                                                dim);
                m++;
                break;
            case 3: 
                medians[first + m] = select_1_3(a[start + startOfGroup],
                                                a[start + startOfGroup + 1],
                                                a[start + startOfGroup + 2],
                                                p,
                                                dim);
                m++;
                break;
            case 4: 
                medians[first + m] = select_1_4(a[start + startOfGroup],
                                                a[start + startOfGroup + 1],
                                                a[start + startOfGroup + 2],
                                                a[start + startOfGroup + 3],
                                                p,
                                                dim);
                m++;
                break;
            default:
                cout << "Error: unhandled case in switch: remainingElements = " << remainingElements << endl;
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
        const T* medianOfMedians =
            medians[ partition(medians, first, m, first + m, (m + 1) >> 1, medians, first + m, p, dim, twoThreads) ];
		     
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
        const long INDEX_CUTOFF = 512;
        if (twoThreads && n > INDEX_CUTOFF) {

            // Yes, more than one thread is available, so calculate the relative index of the middle element.
            long middle = (n + 1) >> 1;

            // Search for the index in the lower half of the array a with a child thread.
            future<void> indexFuture = async(launch::async, [&] {
                for (long i = 0; i < middle; ++i) {
                    if (a[start + i] == medianOfMedians) {
                        swap(a, start + i, start + n - 1);
                        break;
                    }
                }
            });

            // Search for the index in the upper half of the array a with the current thread.
            for (long i = middle; i < n - 1; ++i) {
                if (a[start + i] == medianOfMedians) {
                    swap(a, start + i, start + n - 1);
                    break;
                }
            }

            // Wait for the child thread to finish execution.
            try {
                indexFuture.get();
            } catch (const exception& e) {
                cout << "caught exception " << e.what() << endl;
            }
        } else
#endif
        {
            // No, only one thread is available to find the index of the median of medians.
            for (long i = 0; i < n - 1; ++i) {
                if (a[start + i] == medianOfMedians) {
                    swap(a, start + i, start + n - 1);
                    break;
                }
            }
        }
	    	
        // Partition the array a relative to the median of medians into < and > subsets.
        long i = 0;

#ifdef BIDIRECTIONAL_PARTITION
        // Search from both ends of the array a in order to minimize the use of the swap
        // function at the expense of greater use of the superKeyCompare function.
        //
        // NOTE, however, that searching from both ends appears to degrade performance.
        long j = n - 2;
        while (j > i) {
            if ( superKeyCompare(a[start + i], medianOfMedians, p, dim) < 0 ) {
                ++i;
            } else if ( superKeyCompare(a[start + j], medianOfMedians, p, dim) > 0 ) {
                j--;
            } else {
                swap(a, start + i, start + j);
                ++i;
                j--;
            }
        }

        // Ensure that all elements of the < subset are located below a[start + i].
        for ( ; i < n - 1; ++i) {
            if ( superKeyCompare(a[start + i], medianOfMedians, p, dim) > 0 ) {
                break;
            }
        }
#else
        // Search upward from the beginning of the array a in order to minimize the use of
        // the superKeyCompare function at the expense of greater use of the swap function.
        for (long j = 0; j < n - 1; ++j) {
            if ( superKeyCompare(a[start + j], medianOfMedians, p, dim) < 0 ) {
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
		    	
        } else if (k > i + 1) {
		    	
            // The median of medians occupies a position above i, so partition
            // the array elements of the > subset; for this subset, the
            // original kth element is not the kth element of this subset
            // because i + 1 elements are in the < subset.
            return partition(a, start + i + 1, n - i - 1, size, k - i - 1,
                             medians, first, p, dim, twoThreads);
		    	
        } else {
		    	
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
     * calling parameters:
     *
     * reference - a T** of references to each of the (x, y, z, w...) tuples
     * temporary - a scratch array into which to copy references
     * permutation - avector<long> that indications permutation of the partition coordinate
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
    static KdNode *buildKdTree(T** reference, T** temporary,
                               vector<long>& permutation, vector<KdNode*>& kdNodes,
                               const long start, const long end, const long size,
                               const long dim, const long maximumSubmitDepth,
                               const long depth) {

        KdNode* node = NULL;

        // The partition permutes as x, y, z, w... and specifies the most significant key.
        long p = permutation.at(depth);

        if (end == start) {

            // Only one reference was passed to this method, so store it at this level of the tree.
            node = getKdNode(reference, kdNodes, start);

        } else if (end == start + 1) {
				
            // Two references were passed to this method in unsorted order, so store the
            // start reference at this level of the tree and determine whether to store the
            // end reference as the < child or the > child.
            node = getKdNode(reference, kdNodes, start);
            if (superKeyCompare(reference[start], reference[end], p, dim) > 0) {
                node->ltChild = getKdNode(reference, kdNodes, end);
            } else {
                node->gtChild = getKdNode(reference, kdNodes, end);
            }
				
        } else if (end == start + 2) {
				
            // Three references were passed to this method in unsorted order, so compare
            // the three references to determine which reference is the median reference.
            // Store the median reference at this level of the tree, store the smallest
            // reference as the < child and store the largest reference as the > child.
            long mid = start + 1;
            if (superKeyCompare(reference[start], reference[mid], p, dim) < 0) {
                // reference[start] < reference[mid]
                if (superKeyCompare(reference[mid], reference[end], p, dim) < 0) {
                    // reference[start] < reference[mid] < reference[end]
                    node = getKdNode(reference, kdNodes, mid);
                    node->ltChild = getKdNode(reference, kdNodes, start);
                    node->gtChild = getKdNode(reference, kdNodes, end);
                } else {
                    // reference[start] < reference[mid]; reference[end] < reference[mid]
                    if (superKeyCompare(reference[start], reference[end], p, dim) < 0) {
                        // reference[start] < reference[end] < reference[mid]
                        node = getKdNode(reference, kdNodes, end);
                        node->ltChild = getKdNode(reference, kdNodes, start);
                        node->gtChild = getKdNode(reference, kdNodes, mid);
                    } else {
                        // reference[end] < reference[start] < reference[mid]
                        node = getKdNode(reference, kdNodes, start);
                        node->ltChild = getKdNode(reference, kdNodes, end);
                        node->gtChild = getKdNode(reference, kdNodes, mid);
                    }
                }
            } else {
                // reference[mid] < reference[start]
                if (superKeyCompare(reference[start], reference[end], p, dim) < 0) {
                    // reference[mid] < reference[start] < reference[end]
                    node = getKdNode(reference, kdNodes, start);
                    node->ltChild = getKdNode(reference, kdNodes, mid);
                    node->gtChild = getKdNode(reference, kdNodes, end);
                } else {
                    // reference[mid] < reference[start]; reference[end] < reference[start]
                    if (superKeyCompare(reference[mid], reference[end], p, dim) < 0) {
                        // reference[mid] < reference[end] < reference[start]
                        node = getKdNode(reference, kdNodes, end);
                        node->ltChild = getKdNode(reference, kdNodes, mid);
                        node->gtChild = getKdNode(reference, kdNodes, start);
                    } else { 
                        // reference[end] < reference[mid] < reference[start]
                        node = getKdNode(reference, kdNodes, mid);
                        node->ltChild = getKdNode(reference, kdNodes, end);
                        node->gtChild = getKdNode(reference, kdNodes, start);
                    }
                }
            }
				
        } else if (end > start + 2) {
				
            // Four or more references were passed to this method, so calculate the offset
            // of the median element. partition the reference array about its median element,
            // which is the kth element as calculated below.
            long n = end - start + 1;
            long k = (n + 1) >> 1;

            // Build the < branch of the tree with a child thread at as many levels of the
            // tree as possible.  Create the child thread as high in the tree as possible.
            // Are child threads available to build both branches of the tree?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

                // No, child threads are not available, so find the median element then
                // partition the reference array about it.  Store the median element
                // from the reference array in a new KdNode.
                long median = partition(reference, start, n, size,  k, temporary, start, p, dim, false);
                node = getKdNode(reference, kdNodes, median);

                // Recursively build the < branch of the tree with the current thread.
                node->ltChild = buildKdTree(reference, temporary, permutation, kdNodes, start,
                                            median-1, size, dim, maximumSubmitDepth, depth+1);

                // Then recursively build the > branch of the tree with the current thread.
                node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median+1,
                                            end, size, dim, maximumSubmitDepth, depth+1);

            } else {

                // Yes, child threads are available, so find the median element then partition
                // the reference array about it.  Store the median element from the reference
                // array in a new KdNode.
                long median = partition(reference, start, n, size,  k, temporary, start, p, dim, true);
                node = getKdNode(reference, kdNodes, median);

                // Recursively build the < branch of the tree with a child thread.
                // The recursive call to buildKdTree must be placed in a lambda
                // expression because buildKdTree is a template not a function.
                future<KdNode<T>*> buildFuture = async(launch::async, [&] {
                    return buildKdTree(reference, temporary, permutation, kdNodes, start,
                                       median-1, size, dim, maximumSubmitDepth, depth+1);
                });
                
                // And simultaneously build the > branch of the tree with the current thread.
                node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median+1,
                                            end, size, dim, maximumSubmitDepth, depth+1);

                // Wait for the child thread to finish execution.
                try {
                    node->ltChild = buildFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }

        } else if (end < start) {

            // This is an illegal condition that should never occur, so test for it last.
            cout << "error has occurred at depth = " << depth << " : end = " << end
                 << "  <  start = " << start << endl;
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
     * calling parameters:
     *
     * reference - a T** of references to each of the (x, y, z, w...) tuples
     * temporary - a scratch array into which to copy references
     * permutation - avector<long> that indications permutation of the partition coordinate
     * start - start element of the reference array
     * end - end element of the reference array
     * size - the size of the reference array
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     *
     * returns: a KdNode pointer to the root of the k-d tree
     */
private:
    static KdNode *buildKdTreePresorted(T** reference, T** temporary,
                                        vector<long>& permutation, vector<KdNode*>& kdNodes,
                                        const long start, const long end, const long size,
                                        const long dim, const long maximumSubmitDepth) {

        KdNode *node;

        // It is assumed that the reference array has been pre-sorted using the x:y:z super key.
        long depth = 0;
			
        if (end == start) {

            // Only one reference was passed to this method, so store it at this level of the tree.
            node = getKdNode(reference, kdNodes, start);

        } else if (end == start + 1) {
				
            // Two references were passed to this method in sorted order, so store the start
            // element at this level of the tree and store the end element as the > child. 
            node = getKdNode(reference, kdNodes, start);
            node->gtChild = getKdNode(reference, kdNodes, end);
				
        } else if (end == start + 2) {
				
            // Three references were passed to this method in sorted order, so
            // store the median element at this level of the tree, store the start
            // element as the < child and store the end element as the > child.
            node = getKdNode(reference, kdNodes, start + 1);
            node->ltChild = getKdNode(reference, kdNodes, start);
            node->gtChild = getKdNode(reference, kdNodes, end);
				
        } else if (end > start + 2) {
				
            // Four or more references were passed to this method, so use the median element of
            // the pre-sorted reference array to partition the reference array.
            long n = end - start + 1;
            long median = (n + 1) >> 1;
            node = getKdNode(reference, kdNodes, median);

            // Build the < branch of the tree with a child thread at as many levels of the
            // tree as possible.  Create the child thread as high in the tree as possible.
            // Are child threads available to build both branches of the tree?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

                // No, child threads are not available, so recursively build the < branch
                // of the tree with the current thread.
                node->ltChild = buildKdTree(reference, temporary, permutation, kdNodes, start,
                                            median-1, size, dim, maximumSubmitDepth, depth+1);

                // Then recursively build the > branch of the tree with the current thread.
                node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median+1,
                                            end, size, dim, maximumSubmitDepth, depth+1);

            } else {

                // Yes, child threads are available, so recursively build the < branch
                // of the tree with a child thread. The recursive call to buildKdTree
                // must be placed in a lambda expression because buildKdTree is a template
                // not a function.
                future<KdNode<T>*> buildFuture = async(launch::async, [&] {
                    return buildKdTree(reference, temporary, permutation, kdNodes, start,
                                       median-1, size, dim, maximumSubmitDepth, depth+1);
                });

                // And simultaneously build the > branch of the tree with the current thread.
                node->gtChild = buildKdTree(reference, temporary, permutation, kdNodes, median+1,
                                            end, size, dim, maximumSubmitDepth, depth+1);

                // Wait for the child thread to finish execution.
                try {
                    node->ltChild = buildFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }

        } else if (end < start) {

            // This is an illegal condition that should never occur, so test for it last.
            cout << "error has occurred at depth = " << depth << " : end = " << end
                 << "  <  start = " << start << endl;
            exit(1);

        }

        // Return the pointer to the root of the k-d tree.
        return node;
    }

    /*
     * The verifyKdTree function walks the k-d tree and checks that the
     * children of a node are in the correct branch of that node.
     *
     * calling parameters:
     *
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the depth in the k-d tree
     *
     * returns: a count of the number of kdNodes in the k-d tree
     */
private:
    long verifyKdTree(long const dim, long const maximumSubmitDepth, long const depth) const {
            
        long count = 1 ;
        if (tuple == NULL) {
            cout << "point is null!" << endl;
            exit(1);
        }
            
        // The partition cycles as x, y, z, w...
        long p = depth % dim;
            
        if (ltChild != NULL) {
            if (ltChild->tuple[p] > tuple[p]) {
                cout << "child is > node!" << endl;
                exit(1);
            }
            if (superKeyCompare(ltChild->tuple, tuple, p, dim) >= 0) {
                cout << "child is >= node!" << endl;
                exit(1);
            }
        }
        if (gtChild != NULL) {
            if (gtChild->tuple[p] < tuple[p]) {
                cout << "child is < node!" << endl;
                exit(1);
            }
            if (superKeyCompare(gtChild->tuple, tuple, p, dim) <= 0) {
                cout << "child is <= node" << endl;
                exit(1);
            }
        }
            
        // Verify the < branch with a child thread at as many levels of the tree as possible.
        // Create the child thread as high in the tree as possible for greater utilization.
            
        // Is a child thread available to build the < branch?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                
            // No, so verify the < branch with the current thread.
            if (ltChild != NULL) {
                count += ltChild->verifyKdTree(dim, maximumSubmitDepth, depth + 1);
            }
                
            // Then verify the > branch with the current thread.
            if (gtChild != NULL) {
                count += gtChild->verifyKdTree(dim, maximumSubmitDepth, depth + 1);
            }
        } else {
                
            // Yes, so verify the < branch with a child thread. Note that a
            // lambda is required to instantiate the verifyKdTree template.
            future<long> verifyFuture;
            if (ltChild != NULL) {
                verifyFuture = async(launch::async, [&] {
                    return ltChild->verifyKdTree(dim, maximumSubmitDepth, depth + 1);
                });
            }
                
            // And simultaneously verify the > branch with the current thread.
            long gtCount = 0;
            if (gtChild != NULL) {
                gtCount = gtChild->verifyKdTree(dim, maximumSubmitDepth, depth + 1);
            }
                                                      
            // Wait for the child thread to finish execution.
            long ltCount = 0;
            if (ltChild != NULL) {
                try {
                    ltCount = verifyFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
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
     * calling parameters:
     *
     * coordinates - a vector<T*> of references to each of the (x, y, z, w...) tuples
     * numDimensions - the number of dimensions
     * numThreads - the number of threads
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     *
     * returns: a KdNode pointer to the root of the k-d tree
     */
public:
    static KdNode *createKdTree(vector<T*>& coordinates, const long numDimensions,
                                const long numThreads, const long maximumSubmitDepth) {
        
        struct timespec startTime, endTime;

        // Initialize the reference array.  Create a temporary array for use in
        // sorting the references and building the k-d tree.
        startTime = getTime();
        T** reference = new T*[coordinates.size()];
        T** temporary = new T*[coordinates.size()];
        initializeReference( coordinates, reference );
        endTime = getTime();
        double initTime = (endTime.tv_sec - startTime.tv_sec) +
            1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

        // Sort the reference array using multiple threads. Importantly,
        // for compatibility with the 'permutation' vector initialized below,
        // use the first dimension (0) as the leading key of the super key.
        startTime = getTime();
        mergeSortReferenceAscending(reference, temporary, 0, coordinates.size()-1,
                                    0, numDimensions, maximumSubmitDepth, 0);
        endTime = getTime();
        double sortTime = (endTime.tv_sec - startTime.tv_sec) +
            1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

        // Remove references to duplicate coordinates via one pass through the reference array.
        startTime = getTime();
        long end = removeDuplicates( reference, 0, numDimensions, coordinates.size() );
        endTime = getTime();
        double removeTime = (endTime.tv_sec - startTime.tv_sec) +
            1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

        // Start the timer to time building the k-d tree.
        startTime = getTime();

        // Allocate a vector of KdNodes to avoid contention between multiple threads.
        vector<KdNode*> kdNodes(end + 1);
        for (size_t i = 0; i < kdNodes.size(); ++i) {
            kdNodes[i] = new KdNode();
        }

        // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() ).
        long size = coordinates.size();
        long maxDepth = 1;
        while (size > 0) {
            maxDepth++;
            size >>= 1;
        }
        
        // It is unnecessary to compute the partition coordinate upon each recursive call
        // of the buildKdTree function because that coordinate depends only on the depth of
        // recursion, so it may be pre-computed and stored in the 'permutation' vector.
        // Because the partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
        // (for e.g. 4-dimensional data), the leading key of the super key will be 0 at the
        // first level of the nascent tree, consistent with having sorted the reference array
        // above using 0 as the leading key of the super key.
        vector<long> permutation(maxDepth);
        for (size_t i = 0; i < permutation.size(); ++i) {
            permutation.at(i) = i % numDimensions;
        }
    
        // Build the k-d tree with multiple threads if possible.
        KdNode *root = buildKdTreePresorted(reference, temporary, permutation, kdNodes, 0, end,
                                            coordinates.size(), numDimensions, maximumSubmitDepth);
        endTime = getTime();
        double kdTime = (endTime.tv_sec - startTime.tv_sec) +
            1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));

        // Verify the k-d tree and report the number of kdNodes.
        startTime = getTime();
        long numberOfNodes;
        numberOfNodes = root->verifyKdTree(numDimensions, maximumSubmitDepth, 0);
        endTime = getTime();
        double verifyTime = (endTime.tv_sec - startTime.tv_sec) +
            1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
        cout << "Number of nodes = " << numberOfNodes << endl << endl;

        cout << "totalTime = " << fixed << setprecision(2) << initTime + sortTime + removeTime + kdTime + verifyTime
             << "  initTime = " << initTime << "  sortTime = " << sortTime << "  removeTime = " << removeTime
             << "  kdTime = " << kdTime << "  verifyTime = " << verifyTime << endl << endl;

        // Free the reference and temporary arrays.
        free(reference);
        free(temporary);

        // Return the pointer to the root of the k-d tree.
        return root;
    }

    /*
     * The searchKdTree function searches the k-d tree and find the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
     *
     * calling parameters:
     *
     * query - the query point
     * cut - the cutoff distance
     * dim - the number of dimensions
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the depth in the k-d tree
     *
     * returns: a list that contains the kdNodes that lie within the cutoff distance of the query node
     */
public:
    list< KdNode<T> > searchKdTree(T const* query, long const cut, long const dim,
                                   long const maximumSubmitDepth, long const depth) const {

        // The partition cycles as x, y, z, w...
        long p = depth % dim;
            
        // If the distance from the query node to the KdNode is within the cutoff distance
        // in all k dimensions, add the KdNode to a list.
        list< KdNode<T> > result;
        bool inside = true;
        for (long i = 0; i < dim; ++i) {
            if (abs(query[i] - tuple[i]) > cut) {
                inside = false;
                break;
            }
        }
        if (inside) {
            result.push_back(*this); // The push_back function expects a KdNode for a call by reference.
        }
            
        // Search the < branch with a child thread at as many levels of the tree as possible.
        // Create the child thread as high in the tree as possible for greater utilization.
            
        // Is a child thread available to build the < branch?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {
                
            // No, so search the < branch of the k-d tree with the current thread if the partition
            // coordinate of the query point minus the cutoff distance is <= the partition coordinate
            // of the KdNode.  The < branch must be searched when the cutoff distance equals the
            // partition coordinate because the super key may assign a point to either branch of the
            // tree if the sorting or partition coordinate, which forms the most significant portion
            // of the super key, indicates equality.
            if ( ltChild != NULL && (query[p] - cut) <= tuple[p] ) {
                list<KdNode> ltResult =
                    ltChild->searchKdTree(query, cut, dim, maximumSubmitDepth, depth + 1);
                result.splice(result.end(), ltResult); // Can't substitute searchKdTree(...) for ltResult.
            }
                
            // Then search the > branch of the k-d tree with the current thread if the partition
            // coordinate of the query point plus the cutoff distance is >= the partition coordinate
            // of the KdNode.  The > branch must be searched when the cutoff distance equals the
            // partition coordinate because the super key may assign a point to either branch of the
            // tree if the sorting or partition coordinate, which forms the most significant portion
            // of the super key, indicates equality.
            if ( gtChild != NULL && (query[p] + cut) >= tuple[p] ) {
                list<KdNode> gtResult =
                    gtChild->searchKdTree(query, cut, dim, maximumSubmitDepth, depth + 1);
                result.splice(result.end(), gtResult); // Can't substitute searchKdTree(...) for gtResult.
            }
                
        } else {
                
            // Yes, a child thread is available, so search the < branch with a child thread if the
            // partition coordinate of the query point minus the cutoff distance is <= the partition
            // coordinate of the KdNode.
            future<list<KdNode<T> > > searchFuture;
            if ( ltChild != NULL && (query[p] - cut) <= tuple[p] ) {
                searchFuture = async(launch::async, [&] {
                    return ltChild->searchKdTree(query, cut, dim, maximumSubmitDepth, depth + 1);
                });
            }
                
            // And simultaneously search the < branch with the current thread if the partition coordinate
            // of the query point plus the cutoff distance is >= the partition coordinate of the KdNode.
            list< KdNode<T> > gtResult;
            if ( gtChild != NULL && (query[p] + cut) >= tuple[p] ) {
                gtResult = gtChild->searchKdTree(query, cut, dim, maximumSubmitDepth, depth + 1);
            }
                
            list< KdNode<T> > ltResult;
            if ( ltChild != NULL && (query[p] - cut) <= tuple[p] ) {
                try {
                    ltResult = searchFuture.get();
                } catch (const exception& e) {
                    cout << "caught exception " << e.what() << endl;
                }
            }
            result.splice(result.end(), ltResult);
            result.splice(result.end(), gtResult);
        }
            
        return result;
    }

    /*
     * The printTuple function prints one tuple.
     *
     * calling parameters:
     *
     * tuple - the tuple to print
     * dim - the number of dimensions
     */
public:
    static void printTuple(const T* tuple, const long dim)
        {
            cout << "(" << tuple[0] << ",";
            for (long i=1; i<dim-1; ++i) cout << tuple[i] << ",";
            cout << tuple[dim-1] << ")";
        }

    /*
     * The printKdTree function prints the k-d tree "sideways" with the root at the ltChild.
     *
     * calling parameters:
     *
     * dim - the number of dimensions
     * depth - the depth in the k-d tree 
     */
public:
    void printKdTree(const long dim, const long depth) const
        {
            if (gtChild != NULL) {
                gtChild->printKdTree(dim, depth+1);
            }
            for (long i=0; i<depth; ++i) cout << "       ";
            printTuple(tuple, dim);
            cout << endl;
            if (ltChild != NULL) {
                ltChild->printKdTree(dim, depth+1);
            }
        }
};

/*
 * The randomLongInInterval function creates a random long in the interval [min, max].  See
 * http://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
 *
 * calling parameters:
 *
 * min - the minimum long value desired
 * max - the maximum long value desired
 *
 * returns: a random long
 */
long randomLongInInterval(const long min, const long max) {
    return min + (long) ((((double) rand()) / ((double) RAND_MAX)) * (max - min));
}

/* Create a simple k-d tree and print its topology for inspection. */
int main(int argc, char **argv)
{
    struct timespec startTime, endTime;

    // Set the defaults then parse the input arguments.
    long numPoints = 262144;
    long extraPoints = 100;
    long numDimensions = 3;
    long numThreads = 5;
    long searchDistance = 2000000000;
    long maximumNumberOfNodesToPrint = 5;

    for (long i = 1; i < argc; ++i) {
        if ( 0 == strcmp(argv[i], "-n") || 0 == strcmp(argv[i], "--numPoints") ) {
            numPoints = atol(argv[++i]);
            continue;
        }
        if ( 0 == strcmp(argv[i], "-x") || 0 == strcmp(argv[i], "--extraPoints") ) {
            extraPoints = atol(argv[++i]);
            continue;
        }
        if ( 0 == strcmp(argv[i], "-d") || 0 == strcmp(argv[i], "--numDimensions") ) {
            numDimensions = atol(argv[++i]);
            continue;
        }
        if ( 0 == strcmp(argv[i], "-t") || 0 == strcmp(argv[i], "--numThreads") ) {
            numThreads = atol(argv[++i]);
            continue;
        }
        if ( 0 == strcmp(argv[i], "-s") || 0 == strcmp(argv[i], "--searchDistance") ) {
	    
            continue;
        }
        if ( 0 == strcmp(argv[i], "-p") || 0 == strcmp(argv[i], "--maximumNodesToPrint") ) {
            maximumNumberOfNodesToPrint = atol(argv[++i]);
            continue;
        }
        cout << "illegal command-line argument: " <<  argv[i] << endl;
        exit(1);
    }

    // Declare and initialize the coordinates vector and initialize it with (x,y,z,w) tuples
    // in the half-open interval [0, LONG_MAX] where LONG_MAX is defined in limits.h
    // Create extraPoints-1 duplicate coordinates, where extraPoints <= numPoints,
    // in order to test the removal of duplicate points.
    //
    // Note that the tuples are not vectors in order to avoid copying via assignment statements.
    extraPoints = (extraPoints <= numPoints) ? extraPoints : numPoints;
    vector<long*> coordinates(numPoints + extraPoints - 1);
    for (size_t i = 0; i < coordinates.size(); ++i) {
        coordinates[i] = new long[numDimensions];
    }
    for (long i = 0; i < numPoints; ++i) {
        for (long j = 0; j < numDimensions; ++j) {
            coordinates[i][j] = randomLongInInterval(0, LONG_MAX);
        }
    }
    for (long i = 1; i < extraPoints; ++i) {
        for (long j = 0; j < numDimensions; ++j) {
            coordinates[numPoints - 1 + i][j] = coordinates[numPoints - 1 - i][j];
        }
    }

    // Calculate the number of child threads to be the number of threads minus 1, then
    // calculate the maximum tree depth at which to launch a child thread.  Truncate
    // this depth such that the total number of threads, including the master thread, is
    // an integer power of 2, hence simplifying the launching of child threads by restricting
    // them to only the < branch of the tree for some depth in the tree.
    long n = 0;
    if (numThreads > 0) {
        while (numThreads > 0) {
            n++;
            numThreads >>= 1;
        }
        numThreads = 1 << (n - 1);
    } else {
        numThreads = 0;
    }
    long childThreads = numThreads - 1;
    long maximumSubmitDepth = -1;
    if (numThreads < 2) {
        maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
    } else if (numThreads == 2) {
        maximumSubmitDepth = 0;
    } else {
        maximumSubmitDepth = (long) floor( log( (double) childThreads) / log(2.) );
    }
    cout << endl << "Max number of threads = " << numThreads << "  max submit depth = "
         << maximumSubmitDepth << endl << endl;

    // Create the k-d tree.
    KdNode<long>* root = KdNode<long>::createKdTree(coordinates, numDimensions, numThreads, maximumSubmitDepth);

    // Search the k-d tree for the KdNodes that lie within the cutoff distance of the first tuple.
    long* query = new long[numDimensions];
    for (long i = 0; i < numDimensions; ++i) {
        query[i] = coordinates[0][i];
    }
    startTime = getTime();
    list<KdNode<long> > kdList = root->searchKdTree(query, searchDistance, numDimensions, maximumSubmitDepth, 0);
    endTime = getTime();
    double searchTime = (endTime.tv_sec - startTime.tv_sec) +
        1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
    
    cout << "searchTime = " << fixed << setprecision(2) << searchTime << " seconds" << endl << endl;

    cout << kdList.size() << " nodes within " << searchDistance << " units of ";
    KdNode<long>::printTuple(query, numDimensions);
    cout << " in all dimensions." << endl << endl;
    if (kdList.size() != 0) {
        list<KdNode<long> >::iterator it;
        cout << "List of the first " << maximumNumberOfNodesToPrint << " k-d nodes within a "
             << searchDistance << "-unit search distance follows:" << endl << endl;
        for (it = kdList.begin(); it != kdList.end(); it++) {
            KdNode<long>::printTuple(it->getTuple(), numDimensions);
            cout << endl;
            maximumNumberOfNodesToPrint--;
            if (maximumNumberOfNodesToPrint == 0) {
                break;
            }
    	}
    	cout << endl;
    }	
    return 0;
}
