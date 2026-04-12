/*
 * Copyright (c) 2015, 2019, 2020, 2023, 2025, 2026 Russell A. Brown
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
import java.util.List;
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
    
    private long listNodeCount;
    protected KdNode root, head, tail;
    protected int numDimensions;
    protected ExecutorService executor;
    protected int maxSubmitDepth;

    public KdTree(final int numDimensions,
                  final ExecutorService executor,
                  final int maxSubmitDepth)
    {
                  this.numDimensions = numDimensions;
                  this.executor = executor;
                  this.maxSubmitDepth = maxSubmitDepth;

                  listNodeCount = 0L;
                  root = head = tail = null;
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
     * The {@code remove} method removes a {@code KdNode} from the doubly linked list.
     * 
     * @param node - the {@code KdNode} to prepend
     * @return the length of the updated list
     * </p>
     */
    protected long remove(KdNode node)
    {
        if (node == null) {
            System.out.println("\n\nnull node in KdTree.remove\n");
        }
        if (node.prev == null) {
            // Remove the first node from the list.
            head = node.next;
        } else {
            // Remove any other node from the list.
            node.prev.next = node.next;
        }
        if (node.next == null) {
            // Remove the last node from the list.
            tail = node.prev;
        } else {
            // Remove any other node from the list.
            node.next.prev = node.prev;
        }
        return --listNodeCount;
    }

    /**
     * <p>
     * The {@code add} method prepends appends a {@code KdNode}
     * to the doubly linked list.
     * 
     * @param node - the {@code KdNode} to add
     * @return the length of the updated list
     * </p>
     */
    protected long add(KdNode node)
    {
        if (node == null) {
        System.out.println("\n\nnull node in KdTree.add\n");
        }
        if (head == null) {
            // Add the node to an empty list.
            head = tail = node;
            node.prev = node.next = null;
        } else {
            if(Constants.ENABLE_LIST_PREPEND) {
                // Prepend the node to a non-empty list.
                node.prev = null;
                node.next = head;
                head.prev = node;
                head = node;
            } else {
                // Append the node to a non-empty list.
                node.next = null;
                node.prev = tail;
                tail.next = node;
                tail = node;
            }
        }
        return ++listNodeCount;
     }

    /**
     * <p>
     * The {@code addList} method appends or prepends a doubly linked list from
     * a source {@code KdTree} to the doubly linked list of this {@code KdTree}.
     * 
     * @param tree - the source {@code KdTree}
     * @return the length of the updated list
     * </p>
     */
    protected long addList(final KdTree tree)
    {
        if (head == null && (tree == null || (tree != null && tree.head == null))) {
            // Both lists are empty.
            return 0L;
        } else if (head == null) {
            // This tree's list is empty.
            head = tree.head;
            tail = tree.tail;
            listNodeCount = tree.listSize();
            return listNodeCount;
        } else if (tree == null || (tree != null && tree.head == null)) {
            // The source tree's list is empty.
            return listNodeCount;
        } else {
            // Neither list is empty.
            if (Constants.ENABLE_LIST_PREPEND_ALL) {
                // Prepend the source linked list.
                tree.tail.next = head;
                head.prev = tree.tail;
                tree.tail = tail;
            } else {
                // Append the source linked list.
                tail.next = tree.head;
                tree.head.prev = tail;
                tail = tree.tail;
            }
            return ( listNodeCount + tree.listSize() );
        }
    }

    /**
     * <p>
     * The {@code listSize} method returns the length of the {@code KdNode} list.
     * 
     * @return the length of the list
     * </p>
     */
    protected long listSize()
    {
        if (head == null) {
            return 0L;
        }
        return listNodeCount;
    }

    /**
     * <p>
     * The {@code setList} method sets this {@code KdTree}'s {@code KdNode} list
     * by copying the relevant fields from a source {@code KdTree}.
     * 
     * @param tree - the source {@code KdTree}
     * @return the length of the {@code KdNode} list
     * </p>
     */
    protected long setList(final KdTree tree)
    {
        if (tree == null) {
            return 0L;
        } else {
            head = tree.head;
            tail = tree.tail;
            listNodeCount = tree.listNodeCount;
            return listNodeCount;
        }
    }
    /**
     * <p>
     * The {@code verifyList} method checks the integrity of the {@code KdNode} list.
     * 
     * @return the length of the list
     * </p>
     */
    protected long verifyList()
    {
        if (head == null) {
            return 0L;
        }

        // Count the number of KdNodes in the forward direction.
        long headCount = 0L;
        KdNode node = head;
        while (node != null) {
            ++headCount;
            node = node.next;
        }

        // Count the number of KdNodes in the reverse direction.
       long tailCount = 0L;
        node = tail;
        while (node != null) {
            ++tailCount;
            node = node.prev;
        }

        // Compare the counts.
        if (headCount != listNodeCount) {
            throw new RuntimeException("\n\nhead count = " + headCount +
                                       "  != node count = " + listNodeCount + "\n");
        }

        if (tailCount != listNodeCount) {
            throw new RuntimeException("\n\ntail count = " + tailCount +
                                       "  != node count = " + listNodeCount + "\n");
        }

        return listNodeCount;
    }

    /**
     * <p>
     * The {@code getTreeHeight} method returns the height at the root of the {@code KdTree}.
     * 
     * @return the height of the tree
     * </p>
     */
    protected int getTreeHeight()
    {
        if (root == null) {
            return 0;
        }
        return root.height;
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @param permutation - an array that indicates permutation of the reference arrays
     * @param depth - the depth in the k-d tree
     * @return the number of nodes in the k-d tree
     */
    protected int verifyKdTree(final int[] permutation)
    {
        if (root == null) {
            return 0;
        }
        return root.verifyKdTree(permutation, executor, maxSubmitDepth, 0);
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
    * @return the number of nodes in the k-d tree
     */
    protected int verifyKdTree()
    {
        if (root == null) {
            return 0;
        }
        return root.verifyKdTree(numDimensions, 0, executor, maxSubmitDepth, 0);
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
     * @return the size of the {@link java.util.LinkedList List}{@code <}{@link KdNode}{@code >}
     *         instead of returning void that doesn't appear to be
     *         an acceptable return value for the Callable.call method
     */
    protected List<KdNode> searchKdTree(final long[] queryLower,
                                        final long[] queryUpper,
                                        final ExecutorService executor,
                                        final int maximumSubmitDepth,
                                        final int p,
                                        final int depth,
                                        final boolean enableAll)
    {
        if (root == null) {
            return null;
        }

        if (Constants.ENABLE_LINKED_LIST) {
            return root.searchKdTreeLL(queryLower, queryUpper, executor,
                                       maximumSubmitDepth, p, depth, enableAll);
        } else {
            return root.searchKdTreeAL(queryLower, queryUpper, executor,
                                       maximumSubmitDepth, p, depth, enableAll);
        }
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
    protected List<KdNode> searchKdTree(final long[] queryLower,
                                        final long[] queryUpper,
                                        final ExecutorService executor,
                                        final int maximumSubmitDepth,
                                        final int p,
                                        final int depth,
                                        final boolean[] enable)
    {
        if (root == null) {
            return null;
        }

        if (Constants.ENABLE_LINKED_LIST) {
            return root.searchKdTreeLL(queryLower, queryUpper, executor,
                                       maximumSubmitDepth, p, depth, enable);
        } else {
            return root.searchKdTreeAL(queryLower, queryUpper, executor,
                                       maximumSubmitDepth, p, depth, enable);
        }
    }

    /**
     * <p>
     * The {@code findNearestNeighbors} method finds up to M nearest neighbors to the query array
     * and returns them as a list ordered by increasing distance.
     *
     * Calling parameters:
     *
     * @param query - the query array
     * @param numNeighbors - the number M of nearest neighbors to attempt to find
     * @return a {@link java.util.List List}{@code <}{@link Pair}{@code BigInteger, KdNode>>}
     */
    protected List<Paire> findNearestNeighbors(final long[] query,
                                               final int numNeighbors)
    {
        if (root == null) {
            return null;
        } else {
            return root.findNearestNeighbors(query, numNeighbors);
        }
    }
  
    /**
     * <p>
     * The {@code findNearestNeighbors} method finds up to M nearest neighbors to the query array
     * and returns them as a list ordered by increasing distance.
     *
     * Calling parameters:
     *
     * @param query - the query array
     * @param numNeighbors - the number M of nearest neighbors to attempt to find
     * @param enable - an array that specifies which dimensions to search
     * @return a {@link java.util.List List}{@code <}{@link Pair}{@code BigInteger, KdNode>>}
     */
    protected List<Paire> findNearestNeighbors(final long[] query,
                                               final int numNeighbors,
                                               final boolean[] enable)
    {
        if (root == null) {
            return null;
        } else {
            return root.findNearestNeighbors(query, numNeighbors, enable);
        }
    }
  
    /**
     * <p>
     * The {@code findBruteNeighbors} method finds up to M nearest neighbors to the query array
     * via brute force and returns them as a list ordered by increasing distance.
     *
     * Calling parameters:
     *
     * @param query - the query array
     * @param numNeighbors - the number M of nearest neighbors to attempt to find
     * @return a {@link java.util.List List}{@code <}{@link Pair}{@code BigInteger, KdNode>>}
     */
    protected List<Paire> findBruteNeighbors(final long[] query,
                                             final int numNeighbors)
    {
        if (root == null) {
            return null;
        } else {
            return root.findBruteNeighbors(query, numNeighbors);
        }
    }
  
    /**
     * <p>
     * The {@code findBruteNeighbors} method finds up to M nearest neighbors to the query array
     * via brute force and returns them as a list ordered by increasing distance.
     *
     * Calling parameters:
     *
     * @param query - the query array
     * @param numNeighbors - the number M of nearest neighbors to attempt to find
     * @param enable - an array that specifies which dimensions to search
     * @return a {@link java.util.List List}{@code <}{@link Pair}{@code BigInteger, KdNode>>}
     */
    protected List<Paire> findBruteNeighbors(final long[] query,
                                             final int numNeighbors,
                                             final boolean[] enable)
    {
        if (root == null) {
            return null;
        } else {
            return root.findBruteNeighbors(query, numNeighbors, enable);
        }
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
