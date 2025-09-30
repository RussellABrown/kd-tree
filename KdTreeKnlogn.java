/*
 * Copyright (c) 2015, 2019, 2020, 2023 Russell A. Brown
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

import java.lang.System;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * <p>
 * @author Russell A. Brown
 */

/**
 * <p>
 * The {@code KdTreeKnlogn} class contains static factory methods.
 * </p>
 */
public class KdTreeKnlogn
{
    /**
     * <p>
     * The {@code buildKdTreeKnlogn} method builds a k-d tree by recursively
     * partitioning the reference arrays and adding nodes to the tree. These
     * arrays are permuted cyclically for successive levels of the tree so that
     * sorting uses x, y, z, w, etc. as the most significant portion of the
     * sorting or partitioning key.  The contents of the reference arrays are
     * scrambled by each recursive partitioning.
     * </p>
     *
     * @param references - multiple arrays of KdNode[]
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param start - the first element of the reference array
     * @param end - the last element of the reference array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @returns the root of the k-d tree
     */
    private static KdNode buildKdTreeKnlogn(final KdNode[][] references,
                                            final int[][] permutation,
                                            final int start,
                                            final int end,
                                            final ExecutorService executor,
                                            final int maximumSubmitDepth,
                                            final int depth)
    {
        // This KdNode instance stores the median element.
        final KdNode node;

        // The partition cycles as x, y, z, etc.
        final int p = permutation[depth][permutation[0].length - 1];

        // Get the number of dimensions.
        final int dim = permutation[0].length - 2;

        // Obtain the reference array that corresponds to the most significant key.
        final KdNode[] reference = references[permutation[depth][dim]];

        if (end == start)
        {
            // Only one reference was passed to this method,
            // so store that reference at this level of the tree.
            node = reference[end];
            if (Constants.KD_MAP_DYNAMIC) {
                node.height = 1;
            }

        }
        else if (end == start + 1)
        {
            // Two references were passed to this method in sorted order, so store the start
            // element at this level of the tree and store the end element as the > child. 
            node = reference[start];
            node.gtChild = reference[end];
            if (Constants.KD_MAP_DYNAMIC){
                node.gtChild.height = 1;
                node.height = 2;
            }            
        }
        else if (end == start + 2)
        {
            // Three references were passed to this method in sorted order, so
            // store the median element at this level of the tree, store the start
            // element as the < child and store the end element as the > child.
            node = reference[start + 1];
            node.ltChild = reference[start];
            node.gtChild = reference[end];
            if (Constants.KD_MAP_DYNAMIC){
                node.ltChild.height = node.gtChild.height = 1;
                node.height = 2;
            }
            
        }
        else if (end > start + 2)
        {
            // Four or more references were passed to this method, so the
            // median element of the reference array is chosen as the element
            // about which the other reference arrays will be partitioned
            // Avoid overflow when computing the median.
            int median = start + ((end - start) / 2);
            if (median <= start || median >= end) {
                throw new RuntimeException("error in median calculation at depth = " +
                                           depth + " : start = " + start +
                                           "  median = " + median + "  end = " + end);
            }

            // Store the median element of the reference array in the new KdNode.
            node = reference[median];

            // Build both branches with child threads at as many levels of the tree
            // as possible.  Create the child threads as high in the tree as possible.
            //
            // Is a child thread available to build the < branch, and are there 
            // sufficient KdNode instances to justify spawning a child thread?
            if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth
                || end - start < Constants.KNLOGN_CUTOFF)
            {
                // No, a child thread is not available, so one thread will be used.
                // Initialize startIndex=1 so that the 'for' loop that partitions the
                // reference arrays will partition a number of arrays equal to dim.
                int startIndex = 1;

                // If depth < dim-1, copy references[permut[dim]] to references[permut[0]]
                // where permut is the permutation array for this level of the tree.
                // Sort the two halves of references[permut[0]] with p+1 as the most
                // significant key of the super key. Use as the temporary array
                // references[permut[1]] because that array is not used for partitioning.
                // Partition a number of reference arrays equal to the tree depth because
                // those reference arrays are already sorted.
                if (depth < dim - 1) {
                    startIndex = dim - depth;
                    // Ensure that the partition p cycles as x, y, z, w...
                    final int p1 = (p + 1 < dim) ? p + 1 : 0;
                    // Copy and sort the lower half of references[permut[0]]
                    // with the current thread, excluding the median element.
                    copyAndSortLower(references, reference, permutation, start, median,
                                     p1, executor, maximumSubmitDepth, depth);
                    // Copy and sort the upper half of references[permut[0]]
                    // with the current thread, excluding the median element.
                    copyAndSortUpper(references, reference, permutation, median, end,
                                     p1, executor, maximumSubmitDepth, depth);
                }

                // Partition the reference arrays specified by 'startIndex' in
                // a priori sorted order by comparing super keys.  Store the
                // result from references[permut[i]]] in references[permut[i-1]]
                // where permut is the permutation array for this level of the
                // tree, thus permuting the reference arrays. Skip the element
                // of references[permut[i]] that equals the median element that
                // is stored in the new KdNode.
                final long[] tuple = node.tuple;
                for (int i = startIndex; i < dim; ++i) {
                    // Partiion the lower half of the reference array with the current thread,
                    // excluding the median element.
                    scanAndPartitionLower(references, permutation, node, p, i, start, median, depth);
                    // Partition the upper half of the reference array with the current thread,
                    // excluding the median element.
                    scanAndPartitionUpper(references, permutation, node, p, i, median, end, depth);
                }

                // Recursively build the < branch of the tree with the current thread,
                // excluding the median element.
                node.ltChild = buildKdTreeKnlogn(references, permutation, start, median - 1,
                                                 executor, maximumSubmitDepth, depth + 1);

                // Recursively build the > branch of the tree with the current thread,
                // excluding the median element.
                node.gtChild = buildKdTreeKnlogn(references, permutation, median + 1, end,
                                                 executor, maximumSubmitDepth, depth + 1);
            }
            else
            {
                // Yes, child threads are available, so two threads will be used.
                // Initialize startIndex=1 so that the 'for' loop that partitions the
                // reference arrays will partition a number of arrays equal to dim.
                int startIndex = 1;

                // If depth < dim-1, copy references[permut[dim]] to references[permut[0]]
                // where permut is the permutation array for this level of the tree.
                // Sort the two halves of references[permut[0]] with p+1 as the most
                // significant key of the super key. Use as the temporary array
                // references[permut[1]] because that array is not used for partitioning.
                // Partition a number of reference arrays equal to the tree depth because
                // those reference arrays are already sorted.
                if (depth < dim - 1) {
                    startIndex = dim - depth;
                    // Ensure that the partition p cycles as x, y, z, w...
                    final int p1 = (p + 1 < dim) ? p + 1 : 0;
                    // Copy and sort the lower half of references[permut[0]]
                    // with a child thread, excluding the median element.
                    final Future<Void> future =
                        executor.submit( copyAndSortLowerWithThread(references,
                                                                    reference,
                                                                    permutation,
                                                                    start, median,
                                                                    p1, executor,
                                                                    maximumSubmitDepth,
                                                                    depth) );

                    // And simultaneously copy and sort the upper half of references[permut[0]]
                    // with the current thread, excluding the median element.
                    copyAndSortUpper(references, reference, permutation, median, end,
                                     p1, executor, maximumSubmitDepth, depth);

                    // Wait for the result of copying and sorting the lower half
                    // of the references array with the child thread.
                    try {
                        future.get();
                    } catch (Exception e) {
                        throw new RuntimeException( "copy and sort future exception: " + e.getMessage() );
                    }
                }

                // Partition the reference arrays specified by 'startIndex' in
                // a priori sorted order by comparing super keys.  Store the
                // result from references[permut[i]]] in references[permut[i-1]]
                // where permut is the permutation array for this level of the
                // tree, thus permuting the reference arrays. Skip the element
                // of references[permut[i]] that equals the tuple that is stored
                // in the new KdNode.
                for (int i = startIndex; i < dim; ++i) {

                    // Partition the lower half of the reference array with a child thread,
                    // excluding the median element.
                    final Future<Void> future =
                        executor.submit( scanAndPartitionLowerWithThread(references, permutation,
                                                                         node, p, i,
                                                                         start, median, depth) );

                    // And simultaneously partition the upper half of the reference array
                    // with the current thread, excluding the median element.
                    scanAndPartitionUpper(references, permutation, node, p, i, median, end, depth);

                    // Wait for the result of partitioning the lower half
                    // of the references array with the child thread.
                    try {
                        future.get();
                    } catch (Exception e) {
                        throw new RuntimeException( "partition future exception: " + e.getMessage() );
                    }
                }

                // Recursively build the < branch of the tree with a child thread,
                // excluding the median element.
                final Future<KdNode> future =
                    executor.submit( buildKdTreeKnlognWithThread(references, permutation,
                                                                 start, median - 1, executor,
                                                                 maximumSubmitDepth, depth + 1) );

                // And simultaneously build the > branch of the tree with the current thread,
                // excluding the median element.
                node.gtChild = buildKdTreeKnlogn(references, permutation,
                                                 median + 1, end, executor,
                                                 maximumSubmitDepth, depth + 1);

                // Wait for the result of partitioning the lower half
                // of the references array with the child thread.
                try {
                    node.ltChild = future.get();
                } catch (Exception e) {
                    throw new RuntimeException( "build future exception: " + e.getMessage() );
                }
            }

            if (Constants.KD_MAP_DYNAMIC) {
                // Compute the height at this node as the recursion unwinds.
                node.height = KdTreeDynamic.computeHeight(node);
            }

        } else 	if (end < start) {
            
            // This is an illegal condition that should never occur, so test for it last.
            throw new RuntimeException("end < start");
            
        } else {
            
            // This final else block is added to keep the Java compiler from complaining.
            throw new RuntimeException("unknown configuration of  start and end");
        }
        
        return node;
    }
    
    /**
     * <p>
     * Return a {@link java.util.concurrent.Callable Callable} whose call() method executes the
     * {@link KdTreeKnlogn#buildKdTreeKnlogn buildKdTreeKnlogn} method.
     * </p>
     * 
     * @param references - multiple arrays of KdNode[]
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param start - the first element of the reference array
     * @param end - the last element of the reference array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return a {@link KdNode}
     */
    private static Callable<KdNode> buildKdTreeKnlognWithThread(final KdNode[][] references,
                                                                final int[][] permutation,
                                                                final int start,
                                                                final int end,
                                                                final ExecutorService executor,
                                                                final int maximumSubmitDepth,
                                                                final int depth)
    {
        return new Callable<KdNode>() {
            @Override
            public KdNode call() {
                return buildKdTreeKnlogn(references, permutation,
                                         start, end, executor,
                                         maximumSubmitDepth, depth);
            }
        };
    }

    /**
     * <p>
     * The {@code copyAndSortLower) method copies and then partitions
     * the lower half of a reference array. The median element is excluded.
     * </p>
     *
     * @param references - multiple KdNode[] arrays
     * @param reference - one KdNode[] array
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param start - the first element of the reference array
     * @param median - the median element of the reference array
     * @param p1  - the partition that cycles as x, y, z, etc.
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     */
    private static void copyAndSortLower(final KdNode references[][],
                                         final KdNode reference[],
                                         final int[][] permutation,
                                         final int start,
                                         final int median,
                                         final int p1,
                                         final ExecutorService executor,
                                         final int maximumSubmitDepth,
                                         final int depth)
    {
        KdNode[] dst = references[permutation[depth][0]];
        KdNode[] tmp = references[permutation[depth][1]];
        // Copy and sort elements from start up to but not including median.
        for (int i = start; i < median; ++i) {
            dst[i] = reference[i];
        }
        MergeSort.mergeSortReferenceAscending(dst, tmp, start, median - 1, p1,
                                              executor, maximumSubmitDepth, depth + 1);
    }

    /**
     * <p>
     * Return a {@link java.util.concurrent.Callable Callable} whose call() method
     * executes the {@link KdTreeKnlogn#buildKdTreeKnlogn copyAndSortLower} method.
     * </p>
     *
     * @param references - multiple KdNode[] arrays
     * @param reference - one KdNode[] array
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param start - the first element of the reference array
     * @param median - the median element of the reference array
     * @param p1  - the partition that cycles as x, y, z, etc.
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     */
    private static Callable<Void> copyAndSortLowerWithThread(final KdNode references[][],
                                                             final KdNode reference[],
                                                             final int[][] permutation,
                                                             final int start,
                                                             final int median,
                                                             final int p1,
                                                             final ExecutorService executor,
                                                             final int maximumSubmitDepth,
                                                             final int depth)
    {
        return new Callable<Void>() {
            @Override
            public Void call() {
                copyAndSortLower(references, reference, permutation, start, median,
                                 p1, executor, maximumSubmitDepth, depth);
                return null;
            }
        };
    }
    
    /**
     * <p>
     * The {@code copyAndSortUpper) method copies and then partitions
     * the upper half of a reference array. The median element is excluded.
     * </p>
     *
     * @param references - multiple KdNode[] arrays
     * @param reference - one KdNode[] array
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param median - the median element of the reference array
     * @param end - the last element of the reference array
     * @param p1  - the partition that cycles as x, y, z, etc.
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     */
    private static void copyAndSortUpper(final KdNode references[][],
                                         final KdNode reference[],
                                         final int[][] permutation,
                                         final int median,
                                         final int end,
                                         final int p1,
                                         final ExecutorService executor,
                                         final int maximumSubmitDepth,
                                         final int depth)
    {
        KdNode[] dst = references[permutation[depth][0]];
        KdNode[] tmp = references[permutation[depth][1]];
        // Copy and sort elements from one above median up to and including end.
        for (int i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
        }
        MergeSort.mergeSortReferenceAscending(dst, tmp, median + 1, end, p1,
                                              executor, maximumSubmitDepth, depth + 1);
    }

    /**
     * <p>
     * The {@code scanAndPartitionLower} method scans and partitions
     * the lower half of a reference array.
     * </p>
     *
     * @param references - multiple arrays of KdNode[]
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param node - a {@link KdNode} that contains the coordinate about which to partition
     * @param i - selector for one reference array
     * @param start - the first element of the reference array
     * @param median - the median element of the reference array
     * @param p  - the partition that cycles as x, y, z, etc.
     * @param depth - the depth in the k-d tree
     */
    private static void scanAndPartitionLower(final KdNode references[][],
                                              final int[][] permutation,
                                              final KdNode node,
                                              final int p,
                                              final int i,
                                              final int start,
                                              final int median,
                                              final int depth)
    {
        KdNode src[] = references[permutation[depth][i]];
        KdNode dst[] = references[permutation[depth][i - 1]];
        for (int lower = start - 1, upper = median, j = start; j <= median; ++j) {
            final long compare = MergeSort.superKeyCompare(src[j].tuple, node.tuple, p);
            if (compare < 0L) {
                dst[++lower] = src[j];
            } else if (compare > 0L) {
                dst[++upper] = src[j];
            }
        }
    }
    
    /**
     * <p>
     * Return a {@link java.util.concurrent.Callable Callable} whose call() method executes the
     * {@link KdTreeKnlogn#buildKdTreeKnlogn scanAndPartitionLower} method.
     * </p>
     *
     * @param references - multiple arrays of KdNode[]
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param node - a {@link KdNode} that contains the coordinate about which to partition
     * @param i - selector for one reference array
     * @param start - the first element of the reference array
     * @param median - the median element of the reference array
     * @param p  - the partition that cycles as x, y, z, etc.
     * @param depth - the depth in the k-d tree
     */
    private static Callable<Void> scanAndPartitionLowerWithThread(final KdNode references[][],
                                                                  final int[][] permutation,
                                                                  final KdNode node,
                                                                  final int p,
                                                                  final int i,
                                                                  final int start,
                                                                  final int median,
                                                                  final int depth)
    {
        return new Callable<Void>() {
            @Override
            public Void call() {
                scanAndPartitionLower(references, permutation, node, p, i, start, median, depth);
            return null;
            }
        };
    }
    
    /**
     * <p>
     * Scan and partition the upper half of a reference array.
     * </p>
     *
     * @param references - multiple arrays of KdNode[]
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param node - a {@link KdNode} that contains the coordinate about which to partition
     * @param i - selector for one reference array
     * @param median - the median element of the reference array
     * @param end - the last element of the reference array
     * @param p  - the partition that cycles as x, y, z, etc.
     * @param depth - the depth in the k-d tree
     */
    private static void scanAndPartitionUpper(final KdNode references[][],
                                              final int[][] permutation,
                                              final KdNode node,
                                              final int p,
                                              final int i,
                                              final int median,
                                              final int end,
                                              final int depth)
    {
        KdNode src[] = references[permutation[depth][i]];
        KdNode dst[] = references[permutation[depth][i - 1]];
        for (int lower = median, upper = end + 1, k = end; k > median; --k) {
            final long compare = MergeSort.superKeyCompare(src[k].tuple, node.tuple, p);
            if (compare < 0L) {
                dst[--lower] = src[k];
            } else if (compare > 0L) {
                dst[--upper] = src[k];
            }
        }
    }
    
    /*
    * The swap function swaps two elements in an array.
    *
    * calling parameters:
    *
    * a - the array
    * i - the index of the first element
    * j - the index of the second element
    */
    private static void swap(int a[], int i, int j)
    {
        int t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    
    /**
     * <p>
     * The {@code createKdTreeKnlogn} method builds a k-d tree from a KdNode[]
     * where the coordinates of each point are stored in KdNode.tuple
     * </p>
     *  
     * @param kdNodes - a KdNode[]
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @returns a {@code KdTree} instance
     */
    protected static KdTree createKdTreeKnlogn(final KdNode[] kdNodes,
                                               final ExecutorService executor,
                                               final int maximumSubmitDepth,
                                               final int p)
    {
        // Allocate the references arrays including one additional array.
        final int numPoints = kdNodes.length;
        final int numDimensions = kdNodes[0].tuple.length;
        final KdNode[][] references = new KdNode[numDimensions + 1][numPoints];

        // Don't allocate KdNodes instances for the pth references array
        // (where p is the leading dimension) instead of the first
        // references array to permit KdTreeDynamic.balanceSubtree
        // to build a sub-tree whose root node has a non-zero
        // partition coordinate.
        //
        // Copy references from the KdNode instances of the kdNodes array.
        // These references will be re-ordered by the MergeSort methods.
        for (int i = 0; i < numPoints; ++i) {
            references[p][i] = kdNodes[i];
        }

        // Sort the pth references array using multiple threads, where p
        // is the leading dimension. Importantly, or compatibility with the
        // permutation array initialized below, use the pth dimension as
        // the leading key of the super key. Also, only the pth references
        // array has been populated with KdNode instances.
        MergeSort.mergeSortReferenceAscending(references[p], references[numDimensions],
                                              0, numPoints - 1,
                                              p, executor, maximumSubmitDepth, 0);

        // For a dynamic k-d tree, it is unnecessary to sort the reference array.
        // and remove duplicate coordinates, so merely specify the end index,
        // which is an inlusive index instead of a node count.
        int end = numPoints - 1;            

        // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() ).
        int maxDepth = 1;
        int size = numPoints;
        while (size > 0) {
            ++maxDepth;
            size >>= 1;
        }

        // It is unnecessary to compute the partition coordinate upon each recursive call of
        // the buildKdTree function because that coordinate depends only on the depth of
        // recursion, so it may be pre-computed and stored in the permutation array.
        //
        // first create an array the tracks which reference array contains each coordinate.
        // index 0 is for x, index 1 is for y ... index dim (i.e., numDimensions) is the
        // reference array that MergeSort::mergeSortReferenceAscending used as a temporary array.
        final int[] current = new int[numDimensions + 1];
        for (int i = 0;  i <= numDimensions; ++i) {
            current[i] = i;
        }
        final int[] indices = new int[numDimensions + 2];

        // Create a 2D 'permutation' array from the 'indices' array to specify permutation
        // of the reference arrays and of the partition coordinate.
        final int[][] permutation = new int[maxDepth][numDimensions + 2];

        // Fill the permutation array by calculating the permutation of the indices array
        // and the the partition coordinate of the tuple at each depth in the tree.
        for (int depth = 0; depth < maxDepth; ++depth) {
            int p0 = (depth + p) % numDimensions;
            // The last entry of the indices array contains the partition coordinate.
            indices[numDimensions + 1] = p0;
            // The penultimate entry of the indices array specifies the source reference array.
            indices[numDimensions] = current[p0];
            // The first entry of the indices array specifies the temporary reference array.
            indices[0] = current[numDimensions];
            int k = 1;
            // do the partitioning swaps
            for (int i = 1;  i < numDimensions; ++i) {
                int j = (i + p0) % numDimensions;  // this is the coordinate index following the primary
                indices[k] = current[j];           // write the reference index for that coordinate to the indices array
                swap(current, numDimensions, j);   // this keeps track of where the coordinate will have been have been swapped to.
                ++k;
            }
            // Copy the elements from the indices array to the current row of the permutation array.
            for (int i = 0; i < indices.length; ++i) {
                permutation[depth][i] = indices[i];
            }
        }

        // Build the k-d tree via heirarchical multi-threading if possible.
        final KdTree tree = new KdTree();
        tree.root = buildKdTreeKnlogn(references, permutation, 0, end,
                                      executor, maximumSubmitDepth, 0);
                
        // Return the tree.
        return tree;
    }
    
    /**
     * <p>
     * The {@code createKdTree} method builds a k-d tree from a Pair<Long[], String>[]
     * where the coordinates of each point are stored in Pair.key (a Long[]).
     * </p>
     *  
     * @param coordinates - a Pair<key, value>[]
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param nN, iT, sT, rT, kT, vT - single element arrays for returning values by reference
     * @returns a {@code KdTree} instance
     */
    protected static KdTree createKdTreeKnlogn(final Pair[] coordinates,
                                               final ExecutorService executor,
                                               final int maximumSubmitDepth,
                                               long[] nN,
                                               double[] iT,
                                               double[] sT,
                                               double[] rT,
                                               double[] kT,
                                               double[] vT)
    {
        // Allocate the references arrays including one additional array.
        long initTime = System.currentTimeMillis();
        final int numPoints = coordinates.length;
        final int numDimensions = coordinates[0].getKey().length;
        final KdNode[][] references = new KdNode[numDimensions + 1][numPoints];

        // Allocate KdNodes instances for the first references array.
        // Copy references from the KdNode instances of the kdNodes array.
        // These pointers will be re-ordered by the MergeSort methods.
        for (int i = 0; i < numPoints; ++i) {
            references[0][i] = new KdNode(numDimensions);
            for (int j = 0; j < numDimensions; ++j) {
                references[0][i].tuple[j] = coordinates[i].getKey()[j];
            }
            references[0][i].values.add(coordinates[i].getValue());
        }
        initTime = System.currentTimeMillis() - initTime;

        // Sort the first reference array using the first dimension (0) as the most significant
        // key of the super key.
        long sortTime = System.currentTimeMillis();
        MergeSort.mergeSortReferenceAscending(references[0], references[numDimensions], 0, numPoints - 1,
                                              0, executor, maximumSubmitDepth, 0);
        sortTime = System.currentTimeMillis() - sortTime;

        // Remove references to duplicate tuples via one pass through the sorted reference array.
        // The same most significant key of the super key should be used for sort and de-duping.
        long removeTime = System.currentTimeMillis();
        final int end = MergeSort.removeDuplicates(references[0], 0);
        removeTime = System.currentTimeMillis() - removeTime;

        // Determine the maximum depth of the k-d tree, which is log2( coordinates.size() ).
        int maxDepth = 1;
        int size = numPoints;
        while (size > 0) {
        ++maxDepth;
        size >>= 1;
        }

        // It is unnecessary to compute either the permutation of the reference array or
        // the partition coordinate upon each recursive call of the buildKdTree function
        // because both depend only on the depth of recursion, so they may be pre-computed.
        //
        // Because this array is initialized with 0, 1, 2, 3, 0, 1, 2, 3, etc. (for
        // e.g. 4-dimensional data), the leading key of the super key will be 0 at the
        // first level of the nascent tree, consistent with having sorted the reference
        // array above using 0 as the leading key of the super key.
        //
        // Begin by creating an 'indices' array.
        int[] indices = new int[numDimensions + 2];
        for (int i = 0; i < indices.length - 1; ++i) {
            indices[i] = i;
        }

        // Create a 2D 'permutation' array from the 'indices' array to specify permutation
        // of the reference arrays and of the partition coordinate.
        final int[][] permutation = new int[maxDepth][numDimensions + 2];
        final int[] permutationVerify = new int[maxDepth];

        // Fill the permutation array by calculating the permutation of the indices array
        // and the the partition coordinate of the tuple at each depth in the tree.
        for (int i = 0; i < maxDepth; ++i) {
            // The last entry of the indices array contains the partition coordinate.
            indices[numDimensions + 1] = permutationVerify[i] = i % numDimensions;
            // Swap the first and second to the last elements of the indices array.
            swap(indices, 0, numDimensions);
            // Copy the elements of the indices array to one row of the permutation array.
            for (int j = 0; j < numDimensions + 2; ++j) {
                permutation[i][j] = indices[j];
            }
            // Swap the third and second to the last elements of the indices array.
            swap(indices, numDimensions - 1, numDimensions);
        }

        // Build the k-d tree via heirarchical multi-threading if possible.
        long kdTime = System.currentTimeMillis();
        final KdTree tree = new KdTree();
        tree.root = buildKdTreeKnlogn(references, permutation, 0, end,
                                      executor, maximumSubmitDepth, 0);
        kdTime = System.currentTimeMillis() - kdTime;
        
        // Verify the k-d tree via hierarchical multi-threading if possible and report the number of nodes.
        long verifyTime = System.currentTimeMillis();
        nN[0] = tree.verifyKdTree(permutationVerify, executor, maximumSubmitDepth, 0);
        verifyTime = System.currentTimeMillis() - verifyTime;
       
        iT[0] = (double) initTime / 1000.;
        sT[0] = (double) sortTime / 1000.;
        rT[0] = (double) removeTime / 1000.;
        kT[0] = (double) kdTime / 1000.;
        vT[0] = verifyTime / 1000.;
        
        // Return the tree.
        return tree;
    }
        
} // class KdTreeKnlogn
