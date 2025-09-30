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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.TreeSet;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The {@code KdNode} class stores a point of any number of dimensions
 * as well as references to the "less than" and "greater than" sub-trees.
 * </p>
 */
public class KdNode {
    
    protected long[] tuple;
    protected TreeSet<String> values = null;
    protected int height = 0;
    protected KdNode ltChild = null;
    protected KdNode gtChild = null;

    /**
     * <p>
     * KdNode constructor
     * 
     * @param dim - the number of dimensions of the tuple
     * <p>
     */
    public KdNode(final int numDimensions) {
        tuple = new long[numDimensions];
        values = new TreeSet<String>();
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    protected int verifyKdTree(final int[] permutation,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int depth) {

        if (tuple == null) {
            throw new RuntimeException("point is null");
        }

        // Look up the partition.
        final int p = permutation[depth];

        if (ltChild != null) {
            if (ltChild.tuple[p] > tuple[p]) {
                throw new RuntimeException("node is > partition!");
            }
            if (MergeSort.superKeyCompare(ltChild.tuple, tuple, p) >= 0) {
                throw new RuntimeException("node is >= partition!");
            }
        }
        if (gtChild != null) {
            if (gtChild.tuple[p] < tuple[p]) {
                throw new RuntimeException("node is < partition!");
            }
            if (MergeSort.superKeyCompare(gtChild.tuple, tuple, p) <= 0) {
                throw new RuntimeException("node is <= partition!");
            }
        }
        
        // Count this node.
        int count = 1 ;

        // Search the < branch with a child thread at as many levels of the tree as possible.
        // Create the child thread as high in the tree as possible for greater utilization.

        // Is a child thread available to build the < branch?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

            // No, so search the < branch with the current thread.
            if (ltChild != null) {
                count += ltChild.verifyKdTree(permutation, executor, maximumSubmitDepth, depth + 1);
            }

            // Then search the > branch with the current thread.
            if (gtChild != null) {
                count += gtChild.verifyKdTree(permutation, executor, maximumSubmitDepth, depth + 1);
            }
        } else {

            // Yes, so launch a child thread to search the < branch.
            Future<Integer> future = null;
            if (ltChild != null) {
                future = executor.submit( ltChild.verifyKdTreeWithThread(permutation, executor,
                                                                         maximumSubmitDepth, depth + 1) );
            }
        
            // And simultaneously search the > branch with the current thread.
            if (gtChild != null) {
                count += gtChild.verifyKdTree(permutation, executor, maximumSubmitDepth, depth + 1);
            }

            // If a child thread searched the < branch, get the result.
            if (future != null) {
                try {
                    count += future.get();
                } catch (Exception e) {
                    throw new RuntimeException( "future exception: " + e.getMessage() );
                }
            }
        }

        return count;
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @param dim - the number of dimensions
     * @param q - the leading dimension that permutes cyclically
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    protected int verifyKdTree(final int dim,
                               final int q,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int depth) {

        if (tuple == null) {
            throw new RuntimeException("point is null");
        }

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus opeator for p <= dim.
        final int p = (q < dim) ? q : 0;

        if (ltChild != null) {
            if (ltChild.tuple[p] > tuple[p]) {
                throw new RuntimeException("node is > partition!");
            }
            if (MergeSort.superKeyCompare(ltChild.tuple, tuple, p) >= 0) {
                throw new RuntimeException("node is >= partition!");
            }
        }
        if (gtChild != null) {
            if (gtChild.tuple[p] < tuple[p]) {
                throw new RuntimeException("node is < partition!");
            }
            if (MergeSort.superKeyCompare(gtChild.tuple, tuple, p) <= 0) {
                throw new RuntimeException("node is <= partition!");
            }
        }
        
        // Count this node.
        int count = 1 ;

        // Search the < branch with a child thread at as many levels of the tree as possible.
        // Create the child thread as high in the tree as possible for greater utilization.

        // Is a child thread available to build the < branch?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

            // No, so search the < branch with the current thread.
            if (ltChild != null) {
                count += ltChild.verifyKdTree(dim, p + 1, executor, maximumSubmitDepth, depth + 1);
            }

            // Then search the > branch with the current thread.
            if (gtChild != null) {
                count += gtChild.verifyKdTree(dim, p + 1, executor, maximumSubmitDepth, depth + 1);
            }
        } else {

            // Yes, so launch a child thread to search the < branch.
            Future<Integer> future = null;
            if (ltChild != null) {
                future = executor.submit( ltChild.verifyKdTreeWithThread(dim, p + 1, executor,
                                                                         maximumSubmitDepth, depth + 1) );
            }
        
            // And simultaneously search the > branch with the current thread.
            if (gtChild != null) {
                count += gtChild.verifyKdTree(dim, p + 1, executor, maximumSubmitDepth, depth + 1);
            }

            // If a child thread searched the < branch, get the result.
            if (future != null) {
                try {
                    count += future.get();
                } catch (Exception e) {
                    throw new RuntimeException( "future exception: " + e.getMessage() );
                }
            }
        }

        return count;
    }

    /**
     * <p>
     * The {@code verifyKdTreeWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#verifyKdTree verifyKdTree} method.
     * </p>
     * 
     * @param dim - the number of dimensions
     * @param p - the leading dimension that permutes cyclically
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    private Callable<Integer> verifyKdTreeWithThread(final int dim,
                                                     final int p,
                                                     final ExecutorService executor,
                                                     final int maximumSubmitDepth,
                                                     final int depth) {
        
        return new Callable<Integer>() {
            @Override
                public Integer call() {
                return verifyKdTree(dim, p, executor, maximumSubmitDepth, depth);
            }
        };
    }

    /**
     * <p>
     * The {@code verifyKdTreeWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#verifyKdTree verifyKdTree} method.
     * </p>
     * 
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    private Callable<Integer> verifyKdTreeWithThread(final int[] permutation,
                                                     final ExecutorService executor,
                                                     final int maximumSubmitDepth,
                                                     final int depth) {
        
        return new Callable<Integer>() {
            @Override
                public Integer call() {
                return verifyKdTree(permutation, executor, maximumSubmitDepth, depth);
            }
        };
    }

    /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree to find the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
     * </p>
     *
     * @param result - a {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enableAll - enable or disable all of the k dimensions
     * @return the size of the {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     */
    protected int searchKdTree(final LinkedList<KdNode> result,
                               final long[] queryLower,
                               final long[] queryUpper,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int p,
                               final int depth,
                               final boolean enableAll)
    {
        // Create an enable array of all elements equal to enableAll, and call searchKdTree.
        boolean[] enable = new boolean[queryLower.length];
        Arrays.fill(enable, enableAll);
        return searchKdTree(result, queryLower, queryUpper, executor,
                            maximumSubmitDepth, p, depth, enable);
    }

    /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree to find the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
     * </p>
     *
     * @param result - a {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param q - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enable - enable each of the k dimensions on an individual basis
     * @return the size of the {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     *         instead of returning void that doesn't appear to be
     *         an acceptable return value for the Callable.call method
     */
    protected int searchKdTree(final LinkedList<KdNode> result,
                               final long[] queryLower,
                               final long[] queryUpper,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int q,
                               final int depth,
                               final boolean[] enable)
    {
        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus opeator for p <= dim.
        final int p = (q < queryLower.length) ? q : 0;

        // If the distance from the query node is inside the hyper-rectangle
        // in all k dimensions, prepend the k-d node to the results list.
        boolean inside = true;
        for (int i = 0; i < tuple.length; i++) {
            if (tuple[i] < queryLower[i] || tuple[i] > queryUpper[i]) {
                inside = false;
                break;
            }
        }
        if (inside) {
            result.addFirst(this);
        }

        // Determine whether to search the < and > branches of the k-d tree.
        // If the partition dimension is not enabled, both branches are searched.
        final boolean searchLT = ltChild != null
                                            && (MergeSort.superKeyCompare(tuple, queryLower, p) >= 0
                                            || !enable[p]);
        final boolean searchGT = gtChild != null
                                            && (MergeSort.superKeyCompare(tuple, queryUpper, p) <= 0
                                            || !enable[p]);

        // Do both branches require searching and is a child thread available?
        if (searchLT && searchGT && maximumSubmitDepth >= 0 && depth <= maximumSubmitDepth) {

            // Yes, both branches of the tree require searching and a child thread is available,
            // so prepare to search the < branch with a child thread.
            Future<Integer> future = null;

            // Search the < branch with a child thread.
            LinkedList<KdNode> ltResult = new LinkedList<KdNode>();
            future =
                executor.submit( ltChild.searchKdTreeWithThread(ltResult, queryLower, queryUpper, executor,
                                                                maximumSubmitDepth, p+1, depth+1, enable) );

            // Search the > branch with the master thread?
            LinkedList<KdNode> gtResult = new LinkedList<KdNode>();
            if (searchGT) {
                gtChild.searchKdTree(gtResult, queryLower, queryUpper, executor,
                                     maximumSubmitDepth, p+1, depth+1, enable);
            }

            // Get the result of searching the < branch with the child thread.
            try {
                future.get();
            } catch (Exception e) {
                throw new RuntimeException( "future exception: " + e.getMessage() );
            }

            // Append the results of searching the < and > branches to the result (if any) for this KdNode.
            for (KdNode item : ltResult) {
                result.addFirst(item);
            }
            for (KdNode item : gtResult) {
                result.addFirst(item);
            }

        } else {
        
            // No, both branches do not require searching. Search the < branch with the master thread?
            if (searchLT) {
                LinkedList<KdNode> ltResult = new LinkedList<KdNode>();
                ltChild.searchKdTree(ltResult, queryLower, queryUpper, executor,
                                     maximumSubmitDepth, p+1, depth+1, enable);

                // Append the result of searching the < branch to the result (if any) for this KdNode.
                for (KdNode item : ltResult) {
                    result.addFirst(item);
                }
            }

            // Search the > branch with the master thread?
            if (searchGT) {
                LinkedList<KdNode> gtResult = new LinkedList<KdNode>();
                gtChild.searchKdTree(gtResult, queryLower, queryUpper, executor,
                                     maximumSubmitDepth, p+1, depth+1, enable);

                // Append the result of searching the > branch to the result (if any) for this KdNode.
                for (KdNode item : gtResult) {
                    result.addFirst(item);
                }
            }
        }
        return result.size();
    }

    /**
     * <p>
     * The {@code searchKdTreeWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#searchKdTree searchKdTree} method.
     * </p>
     * 
     * @param result - a {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enable - enable search culling on a per-component basis
     * @return the size of the {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     * that contains the k-d nodes that lie within the cutoff distance of the query node
     */
    public Callable<Integer> searchKdTreeWithThread(final LinkedList<KdNode> result,
                                                    final long[] queryLower,
                                                    final long[] queryUpper,
                                                    final ExecutorService executor,
                                                    final int maximumSubmitDepth,
                                                    final int p,
                                                    final int depth,
                                                    final boolean[] enable) {
        
        return new Callable<Integer>() {
            @Override
                public Integer call() {
                    return searchKdTree(result, queryLower, queryUpper, executor,
                                        maximumSubmitDepth, p, depth, enable);
            }
        };
    }

   /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree to find the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
     * </p>
     *
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enableAll - enable or disable all of the k dimensions
     * @return a {@link java.util.List List}{@code <}{@link KdNode}{@code >}
     * that contains the k-d nodes that lie within the cutoff distance of the query node
     */
    protected ArrayList<KdNode> searchKdTree(final long[] queryLower,
                                             final long[] queryUpper,
                                             final ExecutorService executor,
                                             final int maximumSubmitDepth,
                                             final int p,
                                             final int depth,
                                             final boolean enableAll)
    {
        // Create an enable array of all elements equal to enableAll, and call searchKdTree.
        boolean[] enable = new boolean[queryLower.length];
        Arrays.fill(enable, enableAll);
        return searchKdTree(queryLower, queryUpper, executor,
                            maximumSubmitDepth, p, depth, enable);
    }

    /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree to find the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
     * </p>
     *
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param q - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enable - enable search culling on a per-component basis
     * @return a {@link java.util.List List}{@code <}{@link KdNode}{@code >}
     * that contains the k-d nodes that lie within the cutoff distance of the query node
     */
    protected ArrayList<KdNode> searchKdTree(final long[] queryLower,
                                             final long[] queryUpper,
                                             final ExecutorService executor,
                                             final int maximumSubmitDepth,
                                             final int q,
                                             final int depth,
                                             final boolean[] enable) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus opeator for p <= dim.
        final int p = (q < queryLower.length) ? q : 0;

        // If the distance from the query node is inside the hyper-rectangle
        // in all k dimensions, add the k-d node to a list.
        final ArrayList<KdNode> result = new ArrayList<KdNode>();
        boolean inside = true;
        for (int i = 0; i < tuple.length; i++) {
            if (tuple[i] < queryLower[i] || tuple[i] > queryUpper[i]) {
                inside = false;
                break;
            }
        }
        if (inside) {
            result.add(this);
        }

        // Search the < branch of the k-d tree if the partition coordinate of the query point
        // is greater than the lower side of the query hyper-rectangle for this dimension,
        // or if this component of the query is not enabled for region search.
        Future<ArrayList<KdNode>> future = null;
        if (ltChild != null) {
            if (MergeSort.superKeyCompare(tuple, queryLower, p) >= 0L || !enable[p]) {
            
                // Search the < branch with a child thread at as many levels of the tree as possible.
                // Create the child threads as high in the tree as possible for greater utilization.
                // If maxSubmitDepth == -1, there are no child threads.
                if (maximumSubmitDepth > -1 && depth <= maximumSubmitDepth) {
                    future =
                        executor.submit( ltChild.searchKdTreeWithThread(queryLower, queryUpper, executor,
                                                                        maximumSubmitDepth, p+1, depth+1, enable) );
                } else {
                    result.addAll( ltChild.searchKdTree(queryLower, queryUpper, executor,
                                                        maximumSubmitDepth, p+1, depth+1, enable) );
                }
            }
        }

        // Search the > branch of the k-d tree if the partition coordinate of the query point
        // is less than the upper edge of the query hyper-rectangle for this dimension,
        // or if this component of the query is not enabled for region search.
        if (gtChild != null) {
            if (MergeSort.superKeyCompare(tuple, queryUpper, p) <= 0L || !enable[p]) {
                result.addAll( gtChild.searchKdTree(queryLower, queryUpper, executor,
                                                    maximumSubmitDepth, p+1, depth+1, enable) );
            }
        }
        
        // If a child thread searched the < branch, get the result.
        if (future != null) {
            try {
                result.addAll( future.get() );
            } catch (Exception e) {
                throw new RuntimeException( "future exception: " + e.getMessage() );
            }
        }
        return result;
    }

    /**
     * <p>
     * The {@code searchKdTreeWithThread} method returns a
     * {@link java.util.concurrent.Callable Callable} whose call() method executes the 
     * {@link KdNode#searchKdTree searchKdTree} method.
     * </p>
     * 
     * @param queryLower - the query lower bound array
     * @param queryUpper - the query upper bound array
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     * @param enable - enable search culling on a per-component basis
     * @return a {@link java.util.List List}{@code <}{@link KdNode}{@code >}
     * that contains the k-d nodes that lie within the cutoff distance of the query node
     */
    public Callable<ArrayList<KdNode>> searchKdTreeWithThread(final long[] queryLower,
                                                              final long[] queryUpper,
                                                              final ExecutorService executor,
                                                              final int maximumSubmitDepth,
                                                              final int p,
                                                              final int depth,
                                                              final boolean[] enable) {
        
        return new Callable<ArrayList<KdNode>>() {
            @Override
                public ArrayList<KdNode> call() {
                    return searchKdTree(queryLower, queryUpper, executor,
                                        maximumSubmitDepth, p, depth, enable);
            }
        };
    }

    /**
     * <p>
     * The {@code nearestNeighbor} method is used to search the tree for all possible M
     * nearest geometric neighbors by adding them the the NearestNeighborList.  It
     * excludes branches of the tree where it is guaranteed that all the nodes in that
     * branch are farther way than the current farthest node in the NearestNeighborList.
     * </p>
     *
     * @param nnList - Instance of the NearestNeighborList.
     * @param q - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     */
    protected void nearestNeighbor(final NearestNeighborList nnList,
                                   final int q,
                                   final int depth) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus opeator for p <= dim.
        final int p = (q < nnList.query.length) ? q : 0;

        // Compare the query point to the node point.
        long comparison = MergeSort.superKeyCompare(nnList.query, tuple, p);

        // If query < tuple, descend the < branch to the bottom of the tree before adding a point to the
        // nearestNeighborList, which increases the probability that closer nodes to the query point will get added
        // earlier, thus reducing the likelihood of adding more distant points that get kicked out of the list later.
        if (comparison < 0L) {
            if (ltChild != null) {  // if not at the bottom of the tree yet descend the near branch unconditionally
                ltChild.nearestNeighbor(nnList, p+1, depth+1);
            }
            // Check to see if the current node is closer to the query point than the farthest item in the
            // nearestNeighborList or if this component of the query is not enabled for nearest neighbor search.
            // If so, then add the node and descend the farther branch.
            if (tuple[p] - nnList.query[p] <= nnList.curMaxDist || !nnList.enable[p]) {
                nnList.add(this);  // add the current node to the list
                if (gtChild != null) { // and if not at the bottom, descend the far branch
                    gtChild.nearestNeighbor(nnList, p+1, depth+1);
                }
            }
        }
        // If query > tuple, descend the > branch to the bottom of the tree before adding a point to the
        // nearestNeighborList, which increases the probability that closer nodes to the query point will get added
        // earlier, thus reducing the likelihood of adding more distant points that get kicked out of the list later.
        else if (comparison > 0L) {
            if (gtChild != null) {  // if not at the bottom of the tree yet descend the near branch unconditionally
                gtChild.nearestNeighbor(nnList, p+1, depth+1);
            }
            // Check to see if the current node is closer to the query point than the farthest item in the
            // nearestNeighborList or if this component of the query is not enabled for nearest neighbor search.
            // If so, then add the node and descend the farther branch.
            if (nnList.query[p] - tuple[p] <= nnList.curMaxDist || !nnList.enable[p]) {
                nnList.add(this);
                if (ltChild != null) {
                    ltChild.nearestNeighbor(nnList, p+1, depth+1);
                }
            }
        }
        // query == tuple via comparison of super keys, so the query point is identical to the node point.
        // Add the node then descend both branches because (1) an identical point is by definition the
        // closest that won't get kicked out of the list later, and (2) the probability of finding nearest
        // neighbors is equal for both branches of the tree.
        else {
            nnList.add(this);
            if (ltChild != null) {
                ltChild.nearestNeighbor(nnList, p+1, depth+1);
            }
            if (gtChild != null) {
                gtChild.nearestNeighbor(nnList, p+1, depth+1);
            }
        }
        return;
    }

    /**
     * <p>
     * The {@code printKdTree} method prints the k-d tree "sideways" with the root at the left.
     * </p>
     * 
     * @param depth - the depth in the k-d tree
     */
    protected void printKdTree(final int depth) {
        if (gtChild != null) {
            gtChild.printKdTree(depth+1);
        }
        for (int i = 0; i < depth; i++) {
            System.out.print("        ");
        }
        printTuple(tuple);
        System.out.println();
        if (ltChild != null) {
            ltChild.printKdTree(depth+1);
        }
    }

    /**
     * <p>
     * The {@code printTuple} method prints a tuple.
     * </p>
     * 
     * @param tuple - the tuple
     */
   protected static void printTuple(final long[] tuple) {
        System.out.print("(");
        for (int i = 0; i < tuple.length; i++) {
            System.out.print(tuple[i]);
            if (i < tuple.length - 1) {
                System.out.print(", ");
            }
        }
        System.out.print(")");
    }
} // class KdNode
