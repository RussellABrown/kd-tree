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
 * The {@code KdTreeNlogn} class contains static factory methods.
 * </p>
 */
public class KdTreeNlogn
{
    /**
     * <p>
     * The {@code createKdTreeNlogn} method builds a k-d tree from a KdNode[]
     * where the coordinates of each point are stored in KdNode.tuple
     * </p>
     *  
     * @param kdNodes - a KdNode[]
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @returns a {@code KdTree} instance
     */
    protected static KdTree createKdTreeNlogn(final KdNode[] kdNodes,
                                              final ExecutorService executor,
                                              final int maximumSubmitDepth,
                                              final int p)
    {        
        // Create the reference array and initialize it by
        // copying to it each reference in the kdNodes array.
        final int numPoints = kdNodes.length;
        final int numDimensions = kdNodes[0].tuple.length;
        final KdNode[] reference = new KdNode[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            reference[i] = kdNodes[i];
        }

        // Create the temporary array.
        final KdNode[] temporary = new KdNode[numPoints];
           
        // For a dynamic k-d tree, it is unnecessary to sort the reference array.
        // and remove duplicate coordinates, so merely specify the end index.
        int end = numPoints - 1;
            
        // Determine the maximum depth of the k-d tree, which is log2( numPoints ).
        int size = numPoints;
        int maxDepth = 1;
        while (size > 0) {
            maxDepth++;
            size >>= 1;
        }
            
        // It is unnecessary to compute the partition coordinate upon each recursive call
        // of the buildKdTree function because that coordinate depends only on the depth of
        // recursion, so it may be pre-computed.
        //
        // Add the leading dimension p to the pre-computed partition coordinate (modulo
        // the number of dimensions) to permit KdTreeDynamic::balanceSubtree to build
        // a sub-tree whose root node has a non-zero partition coordinate.
        final int[] permutation = new int[maxDepth];
        for (int i = 0; i < permutation.length; ++i) {
            permutation[i] = (i + p) % numDimensions;
        }
            
        // Build the k-d tree with multiple threads if possible. For a dynamic k-d tree,
        // call the KdNode::buildKdTree function instead of KdNode::buildKdTreePresorted.
        final KdTree tree = new KdTree(numDimensions, executor, maximumSubmitDepth);
        tree.root = buildKdTreeNlogn(reference, temporary, permutation, 0,
                                     end, executor, maximumSubmitDepth, 0);
            
        // Return the tree.
        return tree;
    }
    
    /**
     * <p>
     * The {@code createKdTreeNlogn} method builds a k-d tree from a Pair<Long[], String>[]
     * where the coordinates of each point are stored in Pair.key (a Long[]).
     * </p>
     *  
     * @param coordinates - a Pair<key, value>[]
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param nN, iT, sT, rT, kT, vT - single element arrays for returning values by reference
     * @returns a {@code KdTree} instance
     */
    protected static KdTree createKdTreeNlogn(final Pair[] coordinates,
                                              final ExecutorService executor,
                                              final int maximumSubmitDepth,
                                              long[] nN,
                                              double[] iT,
                                              double[] sT,
                                              double[] rT,
                                              double[] kT,
                                              double[] vT)
    {        
        // Declare and initialize the reference array.
        long initTime = System.currentTimeMillis();
        final int numPoints = coordinates.length;
        final int numDimensions = coordinates[0].getKey().length;
        final KdNode[] reference = new KdNode[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            reference[i] = new KdNode(numDimensions);
            for (int j = 0; j < numDimensions; ++j) {
                reference[i].tuple[j] = coordinates[i].getKey()[j];
            }
            reference[i].values.add(coordinates[i].getValue());
        }
        initTime = System.currentTimeMillis() - initTime;
            
        // Sort the reference array using the first dimension (0) as the most significant
        // key of the super key.
        final KdNode[] temporary = new KdNode[numPoints];
        long sortTime = System.currentTimeMillis();
        MergeSort.mergeSortReferenceAscending(reference, temporary, 0, numPoints - 1,
                                              0, executor, maximumSubmitDepth, 0);
        sortTime = System.currentTimeMillis() - sortTime;

        // Remove references to duplicate coordinates via one pass through the reference array.
        // The same most significant key of the super key should be used for sort and de-duping.
        long removeTime = System.currentTimeMillis();
        final int end = MergeSort.removeDuplicates(reference, 0);
        removeTime = System.currentTimeMillis() - removeTime;
            
        // Determine the maximum depth of the k-d tree, which is log2( numPoints ).
        int size = numPoints;
        int maxDepth = 1;
        while (size > 0) {
            maxDepth++;
            size >>= 1;
        }
            
        // It is unnecessary to compute the partition coordinate upon each recursive call
        // of the buildKdTree function because that coordinate depends only on the depth of
        // recursion, so it may be pre-computed.
        final int[] permutation = new int[maxDepth];
        for (int i = 0; i < permutation.length; ++i) {
            permutation[i] = i % numDimensions;
        }
            
        // Build the k-d tree via hierarchical multi-threading if possible. The
        // buildKdTreePresorted method assumes that the reference array has been
        // pre-sorted using the first dimension (0) as the most significant key
        // of the super key. The permutation array guarantees cycling through the
        // keys so that each key is used as the most significant key of the
        // super key, beginning with key=0 at the first level of tree building.
        long kdTime = System.currentTimeMillis();
        final KdTree tree = new KdTree(numDimensions, executor, maximumSubmitDepth);
        tree.root = buildKdTreePresorted(reference, temporary, permutation, 0,
                                         end, executor, maximumSubmitDepth);
        kdTime = System.currentTimeMillis() - kdTime;
            
        // Verify the k-d tree via hierarchical multi-threading if possible.
        long verifyTime = System.currentTimeMillis();
        nN[0] = tree.verifyKdTree(permutation);
        verifyTime = System.currentTimeMillis() - verifyTime;

        // Return the number of nodes and the execution by reference via arrays.
        iT[0] = (double) initTime / Constants.MILLISECONDS_TO_SECONDS;
        sT[0] = (double) sortTime / Constants.MILLISECONDS_TO_SECONDS;
        rT[0] = (double) removeTime / Constants.MILLISECONDS_TO_SECONDS;
        kT[0] = (double) kdTime / Constants.MILLISECONDS_TO_SECONDS;
        vT[0] = (double) verifyTime / Constants.MILLISECONDS_TO_SECONDS;
            
        // Return the tree.
        return tree;
    }
    
    /**
     * <p>
     * The {@code buildKdTreeNlogn} method builds a k-d tree by recursively partitioning the reference
     * array and adding nodes to the tree.  The super key to be used for partitioning is permuted
     * cyclically for successive levels of the tree in order that sorting use x, y, z, etc. as the
     * most significant portion of the super key.
     * </p>
     *
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] for use in partitioning
     * @param permutation - an array that indicates permutation of the partition coordinate
     * @param start - the first element of the reference array
     * @param end - the last element of the reference array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return a {@link KdNode} that represents the root of the tree or subtree
     */
    private static KdNode buildKdTreeNlogn(final KdNode[] reference,
                                           final KdNode[] temporary,
                                           final int[] permutation,
                                           final int start,
                                           final int end,
                                           final ExecutorService executor,
                                           final int maximumSubmitDepth,
                                           final int depth) {

        // This KdNode instance stores the median element.
        final KdNode node;

        // The partition cycles as x, y, z, etc.
        final int p = permutation[depth];

        if (end == start) {

            // Only one reference was passed to this method, so store it at this level of the tree.
            node = reference[start];
            if (Constants.KD_MAP_DYNAMIC) {
                node.height = 1;
            }

        } else if (end == start + 1) {
                
        // Two references were passed to this method in unsorted order, so store the
        // start reference at this level of the tree and determine whether to store the
        // end reference as the < child or the > child.
        node = reference[start];
        if (MergeSort.superKeyCompare(reference[start].tuple, reference[end].tuple, p) < 0L) {
            node.gtChild = reference[end];
            if (Constants.KD_MAP_DYNAMIC) {
                node.height = 2;
                node.gtChild.height = 1;
            }
        } else {
            node.ltChild = reference[end];
            if (Constants.KD_MAP_DYNAMIC) {
                node.height = 2;
                node.ltChild.height = 1;
            }
        }
                
        } else if (end == start + 2) {
                
        // Three references were passed to this method in unsorted order, so compare
        // the three references to determine which reference is the median reference.
        // Store the median reference at this level of the tree, store the smallest
        // reference as the < child and store the largest reference as the > child.
        int mid = start + 1;
        if (MergeSort.superKeyCompare(reference[start].tuple, reference[mid].tuple, p) < 0L) {
            // reference[start] < reference[mid]
            if (MergeSort.superKeyCompare(reference[mid].tuple, reference[end].tuple, p) < 0L) {
            // reference[start] < reference[mid] < reference[end]
            node = reference[mid];
            node.ltChild = reference[start];
            node.gtChild = reference[end];
            } else {
                // reference[start] < reference[mid]; reference[end] < reference[mid]
                if (MergeSort.superKeyCompare(reference[start].tuple, reference[end].tuple, p) < 0L) {
                    // reference[start] < reference[end] < reference[mid]
                    node = reference[end];
                    node.ltChild = reference[start];
                    node.gtChild = reference[mid];
                } else {
                    // reference[end] < reference[start] < reference[mid]
                    node = reference[start];
                    node.ltChild = reference[end];
                    node.gtChild = reference[mid];
                }
            }
        } else {
            // reference[mid] < reference[start]
            if (MergeSort.superKeyCompare(reference[start].tuple, reference[end].tuple, p) < 0L) {
            // reference[mid] < reference[start] < reference[end]
            node = reference[start];
            node.ltChild = reference[mid];
            node.gtChild = reference[end];
            } else {
                // reference[mid] < reference[start]; reference[end] < reference[start]
                if (MergeSort.superKeyCompare(reference[mid].tuple, reference[end].tuple, p) < 0L) {
                    // reference[mid] < reference[end] < reference[start]
                    node = reference[end];
                    node.ltChild = reference[mid];
                    node.gtChild = reference[start];
                } else { 
                    // reference[end] < reference[mid] < reference[start]
                    node = reference[mid];
                    node.ltChild = reference[end];
                    node.gtChild = reference[start];
                }
            }
        }
        if (Constants.KD_MAP_DYNAMIC) {
            node.height = 2;
            node.ltChild.height = 1;
            node.gtChild.height = 1;
        }
                
        } else if (end > start + 2) {
                
        // Four or more references were passed to this method, so partition the reference
        // array about its median element, which is the kth element as calculated below.
        // Store the median element from the reference array in a new k-d node.
        final int n = end - start + 1;
        final int k = (n + 1) >> 1;
        final int median = partition(reference, start, n, k, temporary, start, p);
        node = reference[median];

        // Build the < branch with a child thread at as many levels of the tree as possible.
        // Create the child threads as high in the tree as possible for greater utilization.
        //
        // Is a child thread available to build the < branch, and are there 
        // sufficient KdNode instances to justify spawning a child thread?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth
            || end - start < Constants.NLOGN_CUTOFF)
        {
            // No, so recursively build the < branch of the tree with the current thread.
            node.ltChild = buildKdTreeNlogn(reference, temporary, permutation, start, median - 1,
                                            executor, maximumSubmitDepth, depth + 1);

            // Then recursively build the > branch of the tree with the current thread.
            node.gtChild = buildKdTreeNlogn(reference, temporary, permutation, median + 1, end,
                                            executor, maximumSubmitDepth, depth + 1);
                    
        }
        else
        {
            // Yes, a child thread is available, so recursively build the < branch with a child thread.
            Future<KdNode> future =
                        executor.submit( buildKdTreeNlognWithThread(reference, temporary, permutation,
                                                                    start, median - 1, executor,
                                                                    maximumSubmitDepth, depth + 1) );
                    
            // And simultaneously, recursively build the > branch of the tree with the current thread.
            node.gtChild = buildKdTreeNlogn(reference, temporary, permutation, median + 1, end,
                                            executor, maximumSubmitDepth, depth + 1);
                    
            // Then get the result of building the < branch with the child thread.
            try {
                node.ltChild = future.get();
            } catch (Exception e) {
                throw new RuntimeException( "future exception: " + e.getMessage() );
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
     * The {@code buildKdTreeWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#buildKdTree buildKdTree} method.
     * </p>
     * 
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] for use in partitioning
     * @param permutation - an array that indicates permutation of the partition coordinate
     * @param start - the first element of the reference array
     * @param end - the last element of the reference array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return a {@link KdNode}
     */
    private static Callable<KdNode> buildKdTreeNlognWithThread(final KdNode[] reference,
                                                               final KdNode[] temporary,
                                                               final int[] permutation,
                                                               final int start,
                                                               final int end,
                                                               final ExecutorService executor,
                                                               final int maximumSubmitDepth,
                                                               final int depth)
    {
        return new Callable<KdNode>() {
        @Override
            public KdNode call() {
                return buildKdTreeNlogn(reference, temporary, permutation, start,
                                        end, executor, maximumSubmitDepth, depth);
            }
        };
    }

    /**
     * <p>
     * The {@code buildKdTreePresorted} method builds a k-d tree by using the median of the
     * pre-sorted reference array to partition that array, then calls the {@link buildKdTree}
     * method to recursively partition the reference array.
     * </p>
     *
     * @param reference - a KdNode[]
     * @param temporary - a scratch KdNode[] for use in partitioning
     * @param permutation - an array that indicates permutation of the partition coordinate
     * @param start - the first element of the reference array
     * @param end - the last element of the reference array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @return a {@link KdNode} that represents the root of the tree or subtree
     */
    private static KdNode buildKdTreePresorted(final KdNode[] reference,
                                               final KdNode[] temporary,
                                               final int[] permutation,
                                               final int start,
                                               final int end,
                                               final ExecutorService executor,
                                               final int maximumSubmitDepth) {

        // This KdNode instance stores the median element.
        final KdNode node;

        // This buildKdTreePresorted() method is used for building only the first level of the tree.
        final int depth = 0;
            
        if (end == start) {

        // Only one reference was passed to this method, so store it at this level of the tree.
        node = reference[start];
        if (Constants.KD_MAP_DYNAMIC) {
            node.height = 1;
        }

        } else if (end == start + 1) {
                
        // Two references were passed to this method in sorted order, so store the start
        // element at this level of the tree and store the end element as the > child. 
        node = reference[start];
        node.gtChild = reference[end];
        if (Constants.KD_MAP_DYNAMIC) {
            node.height = 2;
            node.gtChild.height = 1;
        }
                
        } else if (end == start + 2) {
                
        // Three references were passed to this method in sorted order, so
        // store the median element at this level of the tree, store the start
        // element as the < child and store the end element as the > child.
        node = reference[start + 1];
        node.ltChild = reference[start];
        node.gtChild = reference[end];
        if (Constants.KD_MAP_DYNAMIC) {
            node.height = 2;
            node.ltChild.height = 1;
            node.gtChild.height = 1;
        }
                
        } else if (end > start + 2) {
                
        // Four or more references were passed to this method, so use the median element of
        // the pre-sorted reference array to partition the reference array.
        final int n = end - start + 1;
        final int median = (n + 1) >> 1;
        node = reference[median];

        // Build the < branch with a child thread at as many levels of the tree as possible.
        // Create the child threads as high in the tree as possible for greater utilization.

        // Is a child thread available to build the < branch?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

            // No, so recursively build the < branch of the tree with the current thread.
            node.ltChild = buildKdTreeNlogn(reference, temporary, permutation, start, median - 1,
                                            executor, maximumSubmitDepth, depth + 1);

            // Then recursively build the > branch of the tree with the current thread.
            node.gtChild = buildKdTreeNlogn(reference, temporary, permutation, median + 1, end,
                                            executor, maximumSubmitDepth, depth + 1);
                    
        } else {
                    
            // Yes, a child thread is available, so recursively build the < branch with a child thread.
            final Future<KdNode> future =
                        executor.submit( buildKdTreeNlognWithThread(reference, temporary, permutation,
                                                                    start, median - 1, executor,
                                                                    maximumSubmitDepth, depth + 1) );
                    
            // And simultaneously, recursively build the > branch of the tree with the current thread.
            node.gtChild = buildKdTreeNlogn(reference, temporary, permutation, median + 1, end,
                                            executor, maximumSubmitDepth, depth + 1);
                    
            // Then get the result of building the < branch with the child thread.
            try {
                node.ltChild = future.get();
            } catch (Exception e) {
                throw new RuntimeException( "future exception: " + e.getMessage() );
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
     * The {@code partition} method partitions an array of references to (x,y,z,w...)
     * tuples about its kth element and returns the array index of the kth element.
     * See https://yiqi2.wordpress.com/2013/07/03/median-of-medians-selection-algorithm/
     * that contains a bug that causes a java.lang.ArrayIndexOutOfBoundsException.
     * </p>
     * 
     * @param a - a KdNode[]
     * @param start - the start index for the elements to be considered
     * @param n - the number of elements to consider
     * @param k - the element to find
     * @param medians - a scratch KdNode[] for the medians
     * @param first - the first index for the scratch array
     * @param p - the most significant dimension or the partition coordinate
     * @return the index of the kth element in the array about which the array has been partitioned
     */
    private static int partition(final KdNode[] a,
                                 final int start,
                                 final int n,
                                 final int k,
                                 final KdNode[] medians,
                                 final int first,
                                 final int p)
{
    // The number of elements in a group, which must be 5 due to the select_2_5 method below.
    final int GROUP_SIZE = 5;
        
    if (n <=0 || n > a.length) {
        throw new IllegalArgumentException("n = " + n);
    }
    if (k <= 0 || k > n) {
        throw new IllegalArgumentException("k = " + k);
    }
    if (start + n > a.length) {
        throw new IllegalArgumentException("s = " + start + "  n = " + n + "  length = " + a.length);
    }
        
    // This trivial case terminates recursion.
    if ( n == 1 && k == 1 ) {
        return start;
    }

    // Use insertion sort instead of the median of medians algorithm for a small number of elements,
    // via Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
    // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the a array.
    if (n <= Constants.MEDIAN_OF_MEDIANS_CUTOFF) {
        for (int i = start + 1; i <= start + n - 1; ++i) {
            KdNode tmp = a[i];
            int j;
            for (j = i; j > start && MergeSort.superKeyCompare(a[j-1].tuple, tmp.tuple, p) > 0L; j--) {
                a[j] = a[j-1];
            }
            a[j] = tmp;
        }
        return start + k - 1;
    }
        
    // Otherwise, determine how many medians to find.  Round down to count
    // only groups that comprise fully GROUP_SIZE elements.  Any
    // remaining group of elements that doesn't comprise GROUP_SIZE
    // elements will be processed after the following 'for' loop.
    int m = n / GROUP_SIZE;
        
    // Initialize the index of the start of a group then process each group.
    int startOfGroup = 0;
    for (int i = 0; i < m; i++) {
            
        // Find the median of the group of GROUP_SIZE elements via select_2_5.
        medians[first + i] = select_2_5(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        a[start + startOfGroup + 4],
                                        p);
                
        // Update the index of the beginning of the group of GROUP_SIZE elements.
        startOfGroup += GROUP_SIZE;
    }
        
    // Calculate and check the number of remaining elements.
    int remainingElements = n - startOfGroup;
    if ( remainingElements < 0 || remainingElements >= GROUP_SIZE ) {
        throw new RuntimeException("incorrect group calculation");
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
                                        p);
        m++;
        break;
    case 3: 
        medians[first + m] = select_1_3(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        p);
        m++;
        break;
        case 4: 
        medians[first + m] = select_1_4(a[start + startOfGroup],
                                        a[start + startOfGroup + 1],
                                        a[start + startOfGroup + 2],
                                        a[start + startOfGroup + 3],
                                        p);
        m++;
        break;
    default:
        throw new RuntimeException("unhandled case in switch");
    }
            
    // Select the median of medians for partitioning the elements.  Note that (m + 1) >> 1
    // correctly designates the median element as the "kth" element instead of the address
    // of the median element in the medians array.  The medians array must start with element
    // first + m for the next level of recursion to avoid overwriting the median elements.
    // The medians array will have adequate capacity for all of the medians for all of the
    // recursive calls because the initial call to this partition method from the buildKdTreeNlogn
    // method provides a medians array that is the same size as the reference array.  Each
    // recursive call creates the (1 / GROUP_SIZE) fraction of the medians as the call at the
    // prior level of recursion, so the total requirement for storage of medians is the following
    // fraction of the temporary array for the following values of GROUP_SIZE:
    //
    // for GROUP_SIZE = 3, the fraction is 1/3 + 1/9 + 1/27 + 1/81 + ... < 1/2
    // for GROUP_SIZE = 5, the fraction is 1/5 + 1/25 + 1/125 + 1/625 + ... < 1/4
    // for GROUP_SIZE = 7, the fraction is 1/7 + 1/49 + 1/243 + 1/1701 + ... < 1/8
    //
    // Note: it is possible to allocate the medians array locally to this partition() method
    // instead of providing it via a calling parameter to this method; however, because the
    // mergeSort() method requires a temporary array, that array is re-used as the medians array.
    // But local allocation of the medians array appears to promote marginally faster execution.
    final KdNode medianOfMedians =
        medians[ partition(medians, first, m, (m + 1) >> 1, medians, first + m, p) ];
            
        // Find the address of the median of medians and swap it into a[start + n - 1]
        // so that it is not examined during partitioning of the a array.
        for (int i = 0; i < n - 1; ++i) {
            if (a[start + i] == medianOfMedians) {
                swap(a, start + i, start + n - 1);
                break;
            }
        }
        
        // Partition the a array relative to the median of medians into < and > subsets.
        int i = 0;
        for (int j = 0; j < n - 1; ++j) {
            if (MergeSort.superKeyCompare(a[start + j].tuple, medianOfMedians.tuple, p) < 0L) {
                if (j != i) {
                    swap(a, start + j, start + i);
                }
                ++i;
            }
        }

        // Swap the median of medians into a[start + i] between the < and > subsets.
        swap(a, start + i, start + n - 1);

        // k is 1-based but i is 0-based, so compare k to i + 1 and
        // determine which subset (if any) must be partitioned recursively.
        if (k < i + 1) {
            
            // The median of medians occupies a position below i, so partition
            // the array elements of the < subset; for this subset, the
            // original kth element is still the kth element of this subset.
            return partition(a, start, i, k, medians, first, p);
            
        } else if (k > i + 1) {
            
            // The median of medians occupies a position above i, so partition
            // the array elements of the > subset; for this subset, the
            // original kth element is not the kth element of this subset
            // because i + 1 elements are in the < subset.
            return partition(a, start + i + 1, n - i - 1, k - i - 1,
                                medians, first, p);
            
        } else {
            
            // The median of medians occupies a[start + i] because k == i + 1, so no
            // further partitioning is necessary.  Return start + i as the index of
            // the kth element under the definition that start is the zeroth index.
            return start + i;
        }
    }
        
    /**
     * <p>
     * The {@code swap} method swaps two array elements.
     * </p>
     * 
     * @param a - a KdNode[]
     * @param i - the index of the first element
     * @param j - the index of the second element
     */
    private static void swap(final KdNode[] a, final int i, final int j) {
        final KdNode t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
        
    /**
     * <p>
     * The following {@code select_j_k} methods select the jth of k items.
     * Adapted from Chapter 4, "Linear Orderings", of Alexander Stepanov's
     * and Paul McJones' "Elements of Programming", Addison-Wesley, New York,
     * 2009.
     * </p>
     */
    private static KdNode select_0_2(final KdNode a,
                                     final KdNode b,
                                     final int p)
    {
        if (MergeSort.superKeyCompare(a.tuple, b.tuple, p) < 0L) {
            // a < b
            return a;
            } else {
            // b < a
            return b;
        }
    }

    private static KdNode select_1_2(final KdNode a,
                                     final KdNode b,
                                     final int p)
    {
        if (MergeSort.superKeyCompare(a.tuple, b.tuple, p) < 0L) {
            // a < b
            return b;
        } else {
            // b < a
            return a;
        }
    }
        
    private static KdNode select_1_3_ab(final KdNode a,
                                        final KdNode b,
                                        final KdNode c,
                                        final int p)
    {
        if (MergeSort.superKeyCompare(b.tuple, c.tuple, p) < 0L) {
            // a < b < c
            return b;
        } else {
            // a ? c < b
            return select_1_2(a, c, p);
        }
    }

    private static KdNode select_1_3(final KdNode a,
                                     final KdNode b,
                                     final KdNode c,
                                     final int p)
    {
        if (MergeSort.superKeyCompare(a.tuple, b.tuple, p) < 0L) {
            // a < b
            return select_1_3_ab(a, b, c, p);
        } else {
            // b < a
            return select_1_3_ab(b, a, c, p);
        }
    }

    private static KdNode select_1_4_ab_cd(final KdNode a,
                                           final KdNode b,
                                           final KdNode c,
                                           final KdNode d,
                                           final int p)
    {
        if (MergeSort.superKeyCompare(c.tuple, a.tuple, p) < 0L) {
            // c < a < b && a ? d so c is eliminated and a ? d
            return select_0_2(a, d, p);
        } else {
            // a < b ? c < d so a is eliminated and b ? c
            return select_0_2(b, c, p);
        }
    }
        
    private static KdNode select_1_4_ab(final KdNode a,
                                        final KdNode b,
                                        final KdNode c,
                                        final KdNode d,
                                        final int p)
    {
        if (MergeSort.superKeyCompare(c.tuple, d.tuple, p) < 0L) {
            // a < b && c < d
            return select_1_4_ab_cd(a, b, c, d, p);
        } else {
            // a < b && d < c
            return select_1_4_ab_cd(a, b, d, c, p);
        }
    }
        
    private static KdNode select_1_4(final KdNode a,
                                     final KdNode b,
                                     final KdNode c,
                                     final KdNode d,
                                     final int p)
    {
        if (MergeSort.superKeyCompare(a.tuple, b.tuple, p) < 0L) {
            // a < b
            return select_1_4_ab(a, b, c, d, p);
        } else {
            // b < a
            return select_1_4_ab(b, a, c, d, p);
        }
    }
        
    private static KdNode select_2_5_ab_cd(final KdNode a,
                                           final KdNode b,
                                           final KdNode c,
                                           final KdNode d,
                                           final KdNode e,
                                           final int p)
    {
        if (MergeSort.superKeyCompare(c.tuple, a.tuple, p) < 0L) {
            // c < a < b && c < d ? e so c is eliminated and a < b && d ? e
            return select_1_4_ab(a, b, d, e, p);
        } else {
            // a < b ? c && c < d ? e && b ? e so a is eliminated and c < d && b ? e
            return select_1_4_ab(c, d, b, e, p);
        }
    }
        
    private static KdNode select_2_5_ab(final KdNode a,
                                        final KdNode b,
                                        final KdNode c,
                                        final KdNode d,
                                        final KdNode e,
                                        final int p)
    {
        if (MergeSort.superKeyCompare(c.tuple, d.tuple, p) < 0L) {
        // a < b && c < d
        return select_2_5_ab_cd(a, b, c, d, e, p);
    } else {
        // a < b && d < c
        return select_2_5_ab_cd(a, b, d, c, e, p);
        }
    }
        
    private static KdNode select_2_5(final KdNode a,
                                     final KdNode b,
                                     final KdNode c,
                                     final KdNode d,
                                     final KdNode e,
                                     final int p)
    {
        if (MergeSort.superKeyCompare(a.tuple, b.tuple, p) < 0L) {
            // a < b
            return select_2_5_ab(a, b, c, d, e, p);
        } else {
            // b < a
            return select_2_5_ab(b, a, c, d, e, p);
        }
    }
        
} // class KdTreeNlogn
