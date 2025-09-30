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
import java.util.LinkedList;
import java.util.TreeSet;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The {@code KdTree} class stores the root KdNode of the k-d tree.
 * </p>
 */
public class KdTree {
    
   protected KdNode root;

    protected int getTreeHeight()
    {
        return KdTreeDynamic.getHeight(root);
    }

    /**
     * <p>
     * The {@code createKdTree} method builds a k-d tree from a KdNode[]
     * where the coordinates of each point are stored in KdNode.tuple
     * </p>
     *  
     * @param kdNodes - a KdNode[]
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param p - the leading dimension that permutes cyclically
     * @returns the root of the k-d tree
     */
    protected static KdTree createKdTree(final KdNode[] kdNodes,
                                         final ExecutorService executor,
                                         final int maximumSubmitDepth,
                                         final int p)
    {
        KdTree tree;
        if (Constants.NLOGN) {
            tree = KdTreeNlogn.createKdTreeNlogn(kdNodes, executor,
                                                 maximumSubmitDepth, p);
        } else {
            tree = KdTreeKnlogn.createKdTreeKnlogn(kdNodes, executor,
                                                   maximumSubmitDepth, p);
        }
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
     * @returns the root of the k-d tree
     */
    protected static KdTree createKdTree(final Pair[] coordinates,
                                         final ExecutorService executor,
                                         final int maximumSubmitDepth,
                                         long[] nN,
                                         double[] iT,
                                         double[] sT,
                                         double[] rT,
                                         double[] kT,
                                         double[] vT)
    {
        KdTree tree;
        if (Constants.NLOGN) {
            tree = KdTreeNlogn.createKdTreeNlogn(coordinates, executor,
                                                 maximumSubmitDepth,
                                                 nN, iT, sT, rT, kT, vT);
        } else {
            tree = KdTreeKnlogn.createKdTreeKnlogn(coordinates, executor,
                                                   maximumSubmitDepth,
                                                   nN, iT, sT, rT, kT, vT);
        }
        return tree;
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
                               final int depth)
    {
        if (root == null) {
            return 0;
        } else {
            return root.verifyKdTree(permutation, executor, maximumSubmitDepth, depth);
        }
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @param dim - the number of dimensions
     * @param p - the leading dimension that permutes cyclically
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    protected int verifyKdTree(final int dim,
                               final int p,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int depth)
    {
        if (root == null) {
            return 0;
        } else {
            return root.verifyKdTree(dim, p, executor, maximumSubmitDepth, depth);
        }
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
        return root.searchKdTree(result, queryLower, queryUpper, executor,
                                 maximumSubmitDepth, p, depth, enableAll);
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
     * @param enable - enable search culling on a per-component basis
     * @return the size of the {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     *         instead of returning void that doesn't appear to be
     *         an acceptable return value for the Callable.call method
     */
    protected int searchKdTree(final LinkedList<KdNode> result,
                               final long[] queryLower,
                               final long[] queryUpper,
                               final ExecutorService executor,
                               final int maximumSubmitDepth,
                               final int p,
                               final int depth,
                               final boolean[] enable)
    {
        return root.searchKdTree(result, queryLower, queryUpper, executor,
                                 maximumSubmitDepth, p, depth, enable);
    }

    /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree and finds the KdNodes
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
        return searchKdTree(queryLower, queryUpper, executor,
                            maximumSubmitDepth, p, depth, enableAll);
    }

    /**
     * <p>
     * The {@code searchKdTree} method searches the k-d tree and finds the KdNodes
     * that lie within a cutoff distance from a query node in all k dimensions.
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
    protected ArrayList<KdNode> searchKdTree(final long[] queryLower,
                                             final long[] queryUpper,
                                             final ExecutorService executor,
                                             final int maximumSubmitDepth,
                                             final int p,
                                             final int depth,
                                             final boolean[] enable)
    {
        return root.searchKdTree(queryLower, queryUpper, executor,
                                 maximumSubmitDepth, p, depth, enable);
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
     * @param p - the leading dimension that permutes cyclically
     * @param depth - the depth in the k-d tree
     */
    protected void nearestNeighbor(final NearestNeighborList nnList,
                                   final int p,
                                   final int depth)
    {
        root.nearestNeighbor(nnList, p, depth);
    }

    /**
     * <p>
     * The {@code contains} method searches the tree for a key-value pair.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[], {@code String}>
     * @return {@code true} if the key was found; otherwise, {@code false}
     * </p>
     */
    protected boolean contains(final Pair coordinate)
    {
        if (root != null) {
            return contains(root, coordinate.getKey(), coordinate.getValue(), 0);
        } else {
            return false;
        }
    }
    
    /**
     * <p>
     * The {@code contains} method searches a subtree for a key-value pair.
     * 
     * @param node - the root {@code KdNode} of the subtree
     * @param key - a {@code long}[] tuple
     * @param value - a {@code String}
     * @param q - the leading dimension that is permuted cyclically
     * @return {@code true} if the key was found; otherwise, {@code false}
     * </p>
     */
    private boolean contains(final KdNode node,
                             final long[] key,
                             final String value,
                             final int q)
    {
        KdNode ptr = node;
        int p = q;
        final int dim = node.tuple.length;
        while ( ptr != null ) {
            if (MergeSort.superKeyCompare(key, ptr.tuple, p) < 0L) {
                ptr = ptr.ltChild;
            } else if (MergeSort.superKeyCompare(key, ptr.tuple, p) > 0L) {
                ptr = ptr.gtChild;
            } else {
                // found the key, so check for the value
                return ptr.values.contains(value);
            }
            // Permute the most significant dimension p cyclically using
            // a fast alternative to the modulus operator for p <= dim.
            ++p;
            p = (p < dim) ? p : 0;
        }
        return false; // didn't find the key
    }
    
    /**
     * <p>
     * The {@code printKdTree} method prints the k-d tree "sideways" with the root at the left.
     * </p>
     * 
     * @param depth - the depth in the k-d tree
     */
    protected void printKdTree(final int depth)
    {
        root.printKdTree(depth);
    }

    /**
     * <p>
     * The {@code printTuple} method prints a tuple.
     * </p>
     * 
     * @param p - the tuple
     */
   protected static void printTuple(final long[] tuple)
   {
        KdNode.printTuple(tuple);
   }

} // class KdTree
