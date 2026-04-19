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

/**
 * <p>
 * @author Russell A. Brown
 */

import java.util.concurrent.ExecutorService;
import java.util.TreeSet;

/**
 * <p>
 * The {@code KdTreeDynamic} class extends the (@code KdTree} class and
 * accessess KdTree.root so that an instance of KdTreeDynamic
 * may use KdTree member functions that require access to root.
 * </p>
 */
public class KdTreeDynamic extends KdTree {

    private int nodeCount;
    private boolean inserted, erased, changed;
    protected long[] insertionHistogramDyn, deletionHistogramDyn;
    protected KdNode insertedNode, deletedNode;

    // This constructor does not set KdTree.root
    KdTreeDynamic(final int numDimensions,
                  final ExecutorService executor,
                  final int maxSubmitDepth)
    {
        super(numDimensions, executor, maxSubmitDepth);

        nodeCount = 0;
        inserted = erased = changed = false;
        insertedNode = deletedNode = null;

        // Create histograms to count rebalancing operations;
        // by default, each array element is initialized to 0L.
        insertionHistogramDyn = new long[Constants.MAX_POWER_OF_2];
        deletionHistogramDyn = new long[Constants.MAX_POWER_OF_2];
   }

    // This constructor sets KdTree.root
    KdTreeDynamic(final int numDimensions,
                  final ExecutorService executor,
                  final int maxSubmitDepth,
                  final KdTree tree)
    {
        super(numDimensions, executor, maxSubmitDepth);
        root = tree.root;

        nodeCount = 0;
        inserted = erased = changed = false;
        insertedNode = deletedNode = null;

        // Create histograms to count rebalancing operations.
        if (Constants.ENABLE_HISTOGRAMS) {
            insertionHistogramDyn = new long[Constants.MAX_POWER_OF_2];
            deletionHistogramDyn = new long[Constants.MAX_POWER_OF_2];
        }
    }

    /**
     * <p>
     * The {@code createKdTree} method builds a k-d tree from a list
     * of {@code KdNode}s where the list is provided by a {@code KdTree}
     * and where the coordinates of each point are stored in KdNode.tuple
     * </p>
     *  
     * @param kdTree - a KdTree that holds a list of KdNodes
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @return the root of the k-d tree
     */
    protected static KdTreeDynamic createKdTree(final KdTreeDynamic tree,
                                                final ExecutorService executor,
                                                final int maximumSubmitDepth,
                                                final long[] histogram)
    {
        if (tree == null) {
            throw new RuntimeException("\n\nk-d tree is null in createKdTree\n");
        }
        
        // Create the k-d tree.
        KdTreeDynamic newTree = new KdTreeDynamic(tree.root.tuple.length,
                                                  executor,
                                                  maximumSubmitDepth,
                                                  KdTree.createKdTree(tree,
                                                                      executor,
                                                                      maximumSubmitDepth,
                                                                      histogram));

        // Record the number of k-d nodes in the tree.
        newTree.nodeCount = tree.listSize();

        return newTree;
    }

    /**
     * <p>
     * The {@code createKdTree1to3} method builds a k-d tree from a list
     * of no more than 3 {@code KdNode}s where the list is provided by a
     * {@code KdTree} and where the coordinates of each point are stored
     * in KdNode.tuple
     * </p>
     *  
     * @param kdTree - a KdTree that holds a list of KdNodes
     * @param executor - a {@link java.util.concurrent.ExecutorService ExecutorService}
     * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @return the root of the k-d tree
     */
    protected static KdTreeDynamic createKdTree1to3(final KdTree tree,
                                                    final ExecutorService executor,
                                                    final int maximumSubmitDepth,
                                                    final long[] histogram)
    {
        if (tree == null) {
            throw new RuntimeException("\n\nk-d tree is null in createKdTree1to3\n");
        }
        
        // Create the k-d tree.
        KdTreeDynamic newTree = new KdTreeDynamic(tree.root.tuple.length,
                                                  executor,
                                                  maximumSubmitDepth,
                                                  KdTree.createKdTree1to3(tree,
                                                                          executor,
                                                                          maximumSubmitDepth,
                                                                          histogram));

        // Record the number of k-d nodes in the tree.
        newTree.nodeCount = tree.listSize();

        return newTree;
    }

    /**
     * <p>
     * The {@code treeSize} method returns the size of the {@code KdTreeDynamic}.
     * 
     * @return the tree size if the tree is not null; otherwise, 0
     * </p>
     */
    protected static int getSize(final KdTreeDynamic tree)
    {
        if (tree == null) {
            return 0;
        }
        return tree.nodeCount;
    }

    /**
     * <p>
     * The {@code isEmpty} checks for an empty {@code KdTreeDynamic}.
     * 
     * @return {@code true} if the tree is empty; otherwise {@code false}.
     * </p>
     */
    protected boolean isEmpty()
    {
        return (nodeCount == 0);
    }

    /**
     * <p>
     * The {@code insert} method adds a key-value pair to the tree.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code true} if the key was inserted; otherwise, {@code false}
     * </p>
     */
    protected boolean insert(final Pair coordinate)
    {
        if (coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null in KdTreeDynamic.insert\n");
        }

        inserted = changed = false;
        insertedNode = null;
        if (root != null) {
            root = insert(root, coordinate.getKey(), coordinate.getValue(), insertionHistogramDyn, 0);
            if (insertedNode != null) {
                // A node was inserted into the tree, so count it and add it to the doubly linked list.
                add(insertedNode);
                ++nodeCount;
            }

            // Has the height changed due to an insertion?
            if (changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( !Constants.ENABLE_INSERTION_REBALANCE || isBalanced(root) ) {
                    root.height = computeHeight(root);
                } else {
                    root = rebuildSubTree(root, insertionHistogramDyn, 0);
                }
            }
        } else {
            // Insert the root, count it as an inserted node, add it to the doubly linked list,
            // and increment the histogram element.
            root = insertedNode = new KdNode(coordinate);
            add(insertedNode);
            inserted = changed = true;
            ++nodeCount;
            incrementHistogram(insertionHistogramDyn, 1); // Built a 1-node sub-tree.
        }

        return inserted;
    }

    /**
     * <p>
     * The {@code insert} method adds a key-value pair to a subtree.
     * 
     * @param node - the root {@code KdNode} of the subtree
     * @param key - a {@code long}[] tuple
     * @param value - a {@code String}
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param q - the leading dimension that is permuted cyclically
     * @return the (@code KdNode) root of the subtree
     * </p>
     */
    private KdNode insert(final KdNode node,
                          final long[] key,
                          final String value,
                          final long[] histogram,
                          final int q)
{
        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        final int dim = node.tuple.length;
        final int p = (q < dim) ? q : 0;

        // The return value because the node argument to this function is read only.
        KdNode nodePtr = node;

        // Determine which child to search recursively for an insertion point.
        if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) < 0L) {
            if (nodePtr.ltChild != null) {
                nodePtr.ltChild = insert(nodePtr.ltChild, key, value, histogram, p+1);
            } else {
                // The tree does not contain the key, so create a KdNode instance
                // to insert the pair, and record the reference to that instance.
                nodePtr.ltChild = insertedNode = new KdNode(key, value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) > 0L) {
            if (nodePtr.gtChild != null) {
                nodePtr.gtChild = insert(nodePtr.gtChild, key, value, histogram, p+1);
            } else {
                // The tree does not contain the key, so create a KdNode instance
                // to insert the pair, and record the reference to that instance.
                nodePtr.gtChild = insertedNode = new KdNode(key, value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else {
            // The tree contains the key, so insert the value, even though the
            // values set might already contain it. The tree height doesn't change.
            nodePtr.values.add(value);
            inserted = true;
            changed = false;
        }

        // Has the height changed due to an insertion?
        if (changed) {
            // Yes, the height has changed, so if the subtree
            // rooted at this node remains balanced, compute
            // the height at this node; otherwise rebuild the
            // subtree to rebalance it, which also computes
            // the height. Because rebuilding the subtree
            // recycles its nodes, the node argument to this
            // insert function might no longer specify the
            // root of the subtree.
            if ( !Constants.ENABLE_INSERTION_REBALANCE || isBalanced(nodePtr) ) {
                nodePtr.height = computeHeight(nodePtr);
            } else {
                nodePtr = rebuildSubTree(nodePtr, histogram, p);
            }
        }
        return nodePtr;
    }

    /**
     * <p>
     * The {@code erase} method removes a key-value pair from the tree.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code true} if the key was erased; otherwise, {@code false}
     * </p>
     */
    protected boolean erase(final Pair coordinate)
    {
        if (coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null in KdTreeDynamic.erase\n");
        }

        erased = changed = false;
        deletedNode = null;
        if (root != null) {
            root = erase(root, coordinate.getKey(), coordinate.getValue(), deletionHistogramDyn, 0);
            if (deletedNode != null) {
                // A node was deleted, so count it and remove it from the doubly linked list.
                remove(deletedNode);
                --nodeCount;
            }

            // Has the height changed due to an erasure?
            if (root != null && changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(root) ) {
                    root.height = computeHeight(root);
                } else {
                    root = rebuildSubTree(root, deletionHistogramDyn, 0);
                }
            }
        }
        return erased;
    }

    /**
     * <p>
     * The {@code erase} method removes a key-value pair from a subtree.
     * 
     * @param node - the root {@code KdNode} of the subtree
     * @param key - a {@code long}[] tuple
     * @param value - a {@code String}
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param q - the leading dimension that is permuted cyclically
     * @return the (@code KdNode) root of the subtree
     * </p>
     */
    private KdNode erase(final KdNode node,
                         final long[] key,
                         final String value,
                         final long[] histogram,
                         final int q)
    {
        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        final int dim = node.tuple.length;
        final int p = (q < dim) ? q : 0;

        // The return value because the node argument to this function is read only.
        KdNode nodePtr = node;

        // Determine which child to search recursively for a deletion point.
        if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) < 0L) {
            if (nodePtr.ltChild != null) {
                nodePtr.ltChild = erase(nodePtr.ltChild, key, value, histogram, p+1);

                // Has the height changed due to an erasure? 
                if (erased) {
                    // Yes, the height has changed, so if the subtree
                    // rooted at this node remains balanced, compute
                    // the height at this node; otherwise rebuild the
                    // subtree to rebalance it, which also computes
                    // the height. Because rebuilding the subtree
                    // recycles its nodes, the node argument to this
                    // insert function might no longer specify the
                    // root of the subtree.
                    if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                        nodePtr.height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, histogram, p);
                    }
                }
            } else {
                // The tree does not contain the key.
                erased = changed = false;
            }
        } else if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) > 0L) {
            if (nodePtr.gtChild != null) {
                nodePtr.gtChild = erase(nodePtr.gtChild, key, value, histogram, p+1);
             
                // Has the height changed due to an erasure? 
                if (erased) {
                    // Yes, the height has changed, so if the subtree
                    // rooted at this node remains balanced, compute
                    // the height at this node; otherwise rebuild the
                    // subtree to rebalance it, which also computes
                    // the height. Because rebuilding the subtree
                    // recycles its nodes, the node argument to this
                    // insert function might no longer specify the
                    // root of the subtree.
                    if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                        nodePtr.height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, histogram, p);
                    }
                }
            } else {
                // The tree does not contain the key.
                erased = changed = false;
            }
        } else {
            // The tree contains the key, so if the values set
            // is not empty and the values set does not contain
            // the value, no value is removed and the tree height
            // doesn't change.
            if ( !nodePtr.values.isEmpty() && !nodePtr.values.contains(value) ) {
                erased = changed = false;
            } else {
                // Either the values set is empty or it contains
                // the value. So remove the value from the values set.
                // If the values set is empty, as is the case if the
                // call to this erase method is to erase a replacement
                // node, the call to the TreeSet.remove method is a no-op.
                nodePtr.values.remove(value);

                // Is the values set now empty?
                if ( !nodePtr.values.isEmpty() ) {
                    // The values set is not empty, so a value was erased
                    // but the tree height has not changed.
                    erased = true;
                    changed = false;
                } else {
                    // The values set is empty, so a value was erased and
                    // the tree height will change via removal of a node.
                    erased = changed = true;

                    // Does the node have only a < child?
                    if (nodePtr.ltChild != null && nodePtr.gtChild == null) {
                        // Yes, so the node will be erased from the tree.
                        deletedNode = nodePtr;
                        // Does the subtree rooted at the < child contain <= 3 nodes?
                        //
                        // Checking the height prior to checking the countNodes result
                        // avoids a time-consuming call of countNodes for a large subtree.
                        int nodeCount;
                        if (Constants.ENABLE_1TO3 && nodePtr.ltChild.height <= 3
                            && (nodeCount = countNodes(nodePtr.ltChild)) <= 3)
                        {
                            // Yes, the < child subtree contains <= 3 nodes. So, rebuild
                            // the < child subtree using the leading dimension p at this
                            // level of the tree, and replace the node with its < child.
                            // This approach avoids the need to find a predecessor node.
                            //
                            // Note that rebuildSubTree1to3 computes the height
                            // of the rebuilt < child subtree.
                            nodePtr = rebuildSubTree1to3(nodePtr.ltChild, nodeCount, histogram, p);
                        }
                        else
                        {
                            // No, the < child subtree contains > 3 nodes. So, find
                            // the immediate predecessor node, copy the tuple and the
                            // AvlNode reference from that predecessor node to the
                            // one-child node, swap the values set of that predecessor
                            // node with the values set of the one-child node (whose
                            // values set is currently empty), delete the predecessor
                            // node recursively, and recompute the heights along the
                            // path back to and including the < child.
                            KdNode predecessor = nodePtr.ltChild;
                            predecessor = findPredecessor(nodePtr.ltChild, predecessor, p, p+1);
                            if (Constants.ENABLE_TUPLE_COPY) {
                                System.arraycopy(predecessor.tuple, 0, nodePtr.tuple, 0, dim);
                            } else {
                                nodePtr.tuple = predecessor.tuple;
                            }
                            nodePtr.avlTreeNode = predecessor.avlTreeNode;
                            final TreeSet<String> tmp = predecessor.values;
                            predecessor.values = nodePtr.values;
                            nodePtr.values = tmp;
                            // value is a dummy argument because the predecessor node's
                            // values set is empty.
                            nodePtr.ltChild = erase(nodePtr.ltChild, nodePtr.tuple, value, histogram, p+1);

                            // Assuming that the subtree rooted at the one-child
                            // node was balanced prior to deletion of a node
                            // from the subtree rooted at the < child, the
                            // subtree rooted at the one-child node remains
                            // balanced. The > child subtree is empty, and
                            // deletion of a node from the < child subtree
                            // can decrease but not increase the height of that
                            // subtree. Hence, the balance at the one-child node
                            // can either remain unchanged or decrease.
                            //
                            // Notwithstanding the above reasoning, assume that
                            // the balance may have increased. If the subtree
                            // rooted at the one-child node remains balanced,
                            // compute the height at that node; otherwise,
                            // rebuild the subtree to rebalance it, which also
                            // computes the height. Because rebuilding the subtree
                            // recycles its nodes, the node argument to this
                            // erase function might no longer specify the root
                            // of the subtree.
                            if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                                nodePtr.height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, histogram, p);
                            }
                        }
                    }
                    // Does the node have only a > child?
                    else if (nodePtr.gtChild != null && nodePtr.ltChild == null) {
                        // Yes, so the node will be erased from the tree.
                        deletedNode = nodePtr;
                        // Does the subtree rooted at the > child contain <= 3 nodes?
                        //
                        // Checking the height prior to checking the countNodes result
                        // avoids a time-consuming call of countNodes for a large subtree.
                        int nodeCount;
                        if (Constants.ENABLE_1TO3 && nodePtr.gtChild.height <= 3
                            && (nodeCount = countNodes(nodePtr.gtChild)) <= 3)
                        {
                            // Yes, the > child subtree contains <= 3 nodes. So, rebuild
                            // the > child subtree using the leading dimension p at this
                            // level of the tree, and replace the node with its > child.
                            // This approach avoids the need to find a predecessor node.
                            //
                            // Note that rebuildSubTree1to3 computes the height
                            // of the rebuilt > child subtree.
                            nodePtr = rebuildSubTree1to3(nodePtr.gtChild, nodeCount, histogram, p);
                        }
                        else
                        {
                            // No, the > child subtree contains > 3 nodes. So, find
                            // the immediate successor node, copy the tuple and the
                            // AvlNode reference from that successor node to the
                            // one-child node, swap the values set of that successor
                            // node with the values set of the one-child node (whose
                            // values set is currently empty), delete the successor
                            // node recursively, and recompute the heights along the
                            // path back to and including the > child.
                            KdNode successor = nodePtr.gtChild;
                            successor = findSuccessor(nodePtr.gtChild, successor, p, p+1);
                            if (Constants.ENABLE_TUPLE_COPY) {
                                System.arraycopy(successor.tuple, 0, nodePtr.tuple, 0, dim);
                            } else {
                                nodePtr.tuple = successor.tuple;
                            }
                            nodePtr.avlTreeNode = successor.avlTreeNode;
                            final TreeSet<String> tmp = successor.values;
                            successor.values = nodePtr.values;
                            nodePtr.values = tmp;
                            // value is a dummy argument because the successor node's
                            // values set is empty.
                            nodePtr.gtChild = erase(nodePtr.gtChild, nodePtr.tuple, value, histogram, p+1);

                            // Assuming that the subtree rooted at the one-child
                            // node was balanced prior to deletion of a node
                            // from the subtree rooted at the > child, the
                            // subtree rooted at the one-child node remains
                            // balanced. The < child subtree is empty, and
                            // deletion of a node from the > child subtree
                            // can decrease but not increase the height of that
                            // subtree. Hence, the balance at the one-child node
                            // can either remain unchanged or decrease.
                            //
                            // Notwithstanding the above reasoning, assume that
                            // the balance may have increased. If the subtree
                            // rooted at the one-child node remains balanced,
                            // compute the height at that node; otherwise,
                            // rebuild the subtree to rebalance it, which also
                            // computes the height. Because rebuilding the subtree
                            // recycles its nodes, the node argument to this
                            // erase function might no longer specify the root
                            // of the subtree.
                            if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                                nodePtr.height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, histogram, p);
                            }
                        }
                    }
                    // If the node has no children, set the nodePtr return value
                    // to null so that the node's parent's child reference will
                    // set to null and hence the node will be garbage collected.
                    else if (nodePtr.ltChild == null && nodePtr.gtChild == null) {
                        deletedNode = nodePtr;
                        nodePtr = null;
                    }
                    // The node has two children.
                    else { 
                        // Does the subtree rooted at this two-child node 
                        // have <= 3 nodes, excluding this two-child node?
                        //
                        // Checking the height prior to checking the countNodesSkipRoot result
                        // avoids a time-consuming call of countNodesSkipRoot for a large subtree.
                        int nodeCount;
                        if (Constants.ENABLE_1TO3 && nodePtr.height <= 3
                            && (nodeCount = countNodesSkipRoot(nodePtr)) <= 3)
                        {
                            // Yes, the subtree rooted at this two-child node contains
                            // <= 3 nodes, excluding this two-child node, so rebuild
                            // the subtree excluding this two-child node, thus deleting
                            // this two-child node. This approach avoids the need to
                            // find a predecessor or successor node.
                            deletedNode = nodePtr;
                            nodePtr = rebuildSubTreeSkipRoot1to3(nodePtr, nodeCount, histogram, p);
                        }
                        else
                        {
                            // No, the subtree rooted at this two-child node contains
                            // > 3 nodes, excluding this two-child node, so replace this
                            // two-child node by either its predecessor or its successor.
                            //
                            // If Constants.ENABLE_PREFERRED_NODE is true, select the replacement
                            // node from the taller of the child subtrees.
                            if (Constants.ENABLE_PREFERRED_NODE
                                && (getHeight(nodePtr.ltChild) >= getHeight(nodePtr.gtChild)))
                            {
                                // Find the node with the largest super-key in the
                                // subtree rooted at the < child, which is the
                                // immediate predecessor node. Copy the tuple and the
                                // AvlTree reference from that predecessor node to the
                                // two-child node, swap the values set of that predecessor
                                // node with the values set of the two-child node (whose
                                // values set is currently empty), delete the predecessor
                                // node recursively, and recompute the heights along the
                                // path back to and including the < child.
                                KdNode predecessor = nodePtr.ltChild;
                                predecessor = findPredecessor(nodePtr.ltChild, predecessor, p, p+1);
                                if (Constants.ENABLE_TUPLE_COPY) {
                                    System.arraycopy(predecessor.tuple, 0, nodePtr.tuple, 0, dim);
                                } else {
                                    nodePtr.tuple = predecessor.tuple;
                                }
                                    nodePtr.avlTreeNode = predecessor.avlTreeNode;
                                final TreeSet<String> tmp = predecessor.values;
                                predecessor.values = nodePtr.values;
                                nodePtr.values = tmp;
                                 // value is a dummy argument because the predecessor node's
                                // values set is empty
                                nodePtr.ltChild = erase(nodePtr.ltChild, nodePtr.tuple, value, histogram, p+1);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                                    nodePtr.height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, histogram, p);
                                }
                            }
                            else
                            {
                                // Find the node with the smallest super-key in the
                                // subtree rooted at the > child, which is the
                                // immediate successor node. Copy the tuple and the
                                // AvlTree reference from that successor node to the
                                // two-child node, swap the values set of that successor
                                // node with the values set of the two-child node (whose
                                // values set is currently empty), delete the successor
                                // node recursively, and recompute the heights along the
                                // path back to and including the > child.
                                KdNode successor = nodePtr.gtChild;
                                successor = findSuccessor(nodePtr.gtChild, successor, p, p+1);
                                if (Constants.ENABLE_TUPLE_COPY) {
                                    System.arraycopy(successor.tuple, 0, nodePtr.tuple, 0, dim);
                                } else {
                                    nodePtr.tuple = successor.tuple;
                                }
                                nodePtr.avlTreeNode = successor.avlTreeNode;
                                final TreeSet<String> tmp = successor.values;
                                successor.values = nodePtr.values;
                                nodePtr.values = tmp;
                                // value is a dummy argument because the successor node's
                                // values set is empty.
                                nodePtr.gtChild = erase(nodePtr.gtChild, nodePtr.tuple, value, histogram, p+1);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( !Constants.ENABLE_DELETION_REBALANCE || isBalanced(nodePtr) ) {
                                    nodePtr.height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, histogram, p);
                                }
                            }
                        }
                    }
                }
            }
        }
        return nodePtr;
    }

    /**
     * <p>
     * The {@code findPredecessor} method searches a subtree
     * for the immediate predecessor {@code KdNode}, which is
     * the node with the largest super-key.
     * 
     * @param node - the {@code KdNode} root of the subtree
     * @param predecessor - the current precedecessor {@code KdNode}
     * @param p0 - the leading dimension of the {@code KdNode} to be replaced
     * @param q - the leading dimension that is permuted cyclically
     * @return the updated immediate predecessor {@code KdNode}
     * </p>
     */
    private KdNode findPredecessor(final KdNode node,
                                   final KdNode predecessor,
                                   final int p0,
                                   final int q) {

        // The return value because the predecessor argument to this function is read only.
        KdNode pred = predecessor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        final int dim = node.tuple.length;
        final int p = (q < dim) ? q : 0;

        // Does the leading dimension at this node equal
        // the leading dimension at the node to be replaced?
        if (p == p0) {
            // Yes, the leading dimensions are equal, so if this
            // node has a > child, follow that child; otherwise,
            // this node is a potential predecessor, so check it.
            if (node.gtChild != null) {
                pred = findPredecessor(node.gtChild, pred, p0, p+1);
            } else {
                pred = checkPredecessor(node, pred, p0);
            }
        } else {
            // No, the leading dimensions are not equal; hence, this
            // node is a potential predecessor, so check it, and
            // then inspect both children.
            pred = checkPredecessor(node, pred, p0);
            if (node.ltChild != null) {
                pred = findPredecessor(node.ltChild, pred, p0, p+1);
            }
            if (node.gtChild != null) {
                pred = findPredecessor(node.gtChild, pred, p0, p+1);
            }
        }
        return pred;
    }

    /**
     * <p>
     * The {@code checkPredecessor} checks a potential predecessor
     * (@code KdNode) to see whether it is a larger precedecessor
     * {@code KdNode} than the current predecessor {@code KdNode}.
     * 
     * @param node - a potential prdecessor {@code KdNode}
     * @param predecessor - the current precedecessor {@code KdNode}
     * @param p0 - the leading dimension of the {@code KdNode} to be replaced
     * @return the updated immediate predecessor {@code KdNode}
     * </p>
     */
    private KdNode checkPredecessor(final KdNode node,
                                    final KdNode predecessor,
                                    final int p0) {

        // The return value because the predecessor argument to this function is read only.
        KdNode pred = predecessor;

        // Update the current precedessor if the potential predecessor is larger.
        if (MergeSort.superKeyCompare(node.tuple, pred.tuple, p0) > 0L) {
            pred = node;
        }
        return pred;
    }

    /**
     * <p>
     * The {@code findSuccessor} method searches a subtree
     * for the immediate successor {@code KdNode}, which is
     * the node with the smallest super-key.
     * 
     * @param node - the {@code KdNode} root of the subtree
     * @param successor - the current successor {@code KdNode}
     * @param p0 - the leading dimension of the {@code KdNode} to be replaced
     * @param q - the leading dimension that is permuted cyclically
     * @return the updated immediate successor {@code KdNode}
     * </p>
     */
    private KdNode findSuccessor(final KdNode node,
                                 final KdNode successor,
                                 final int p0,
                                 final int q) {

        // The return value because the successor argument to this function is read only.
        KdNode succ = successor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        final int dim = node.tuple.length;
        final int p = (q < dim) ? q : 0;

        // Does the leading dimension at this node equal
        // the leading dimension at the node to be replaced?
         if (p == p0) {
            // Yes, the leading dimensions are equal, so if this
            // node has a < child, follow that child; otherwise,
            // this node is a potential successor, so check it.
            if (node.ltChild != null) {
                succ = findSuccessor(node.ltChild, succ, p0, p+1);
            } else {
                succ = checkSuccessor(node, succ,  p0);
            }
        } else {
            // No, the leading dimensions are not equal; hence, this
            // node is a potential successor, so check it, and
            // then inspect both children
            succ = checkSuccessor(node, succ, p0);
            if (node.ltChild != null) {
                succ = findSuccessor(node.ltChild, succ, p0, p+1);
            }
            if (node.gtChild != null) {
                succ = findSuccessor(node.gtChild, succ, p0, p+1);
            }
        }
        return succ;
    }

    /**
     * <p>
     * The {@code checkSuccessor} checks a potential successor
     * (@code KdNode) to see whether it is a smaller succcessor
     * {@code KdNode} than the current succcessor {@code KdNode}.
     * 
     * @param node - a potential prdecessor {@code KdNode}
     * @param successor - the current successor {@code KdNode}
     * @param p0 - the leading dimension of the {@code KdNode} to be replaced
     * @return the updated immediate successor {@code KdNode}
     * </p>
     */
    private KdNode checkSuccessor(final KdNode node,
                                  final KdNode successor,
                                  final int p0) {

        // The return value because the successor argument to this function is read only.
        KdNode succ = successor;

        // Update the current successor if the potential successor is smaller.
        if (MergeSort.superKeyCompare(node.tuple, succ.tuple, p0) < 0L) {
            succ = node;
        }
        return succ;
    }

   /**
     * <p>
     * The {@code rebuildSubTree} method rebuilds a subtree
     * rooted at a {@code KdNode} to rebalance it.
     * 
     * @param node - the {@code KdNode} instance
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    protected KdNode rebuildSubTree(final KdNode node,
                                    final long[] histogram,
                                    final int p)
    {
        // Count the nodes in the subtree, allocate an array of KdNode
        // references, and assign a KdNode reference to each array element.
        final int count = countNodes(node);
        final KdNode[] kdNodes = new KdNode[count];
        getSubTree(node, kdNodes);
        
        // If Constants.ENABLE_1TO3 is true, check whether
        // the subtree contains 3 KdNodes or fewer.
        if (Constants.ENABLE_1TO3 && count <= 3) {
            // The subtree contains 3 nodes or fewer,
            // so rebuild the subtree by explicitly
            // comparing the nodes' super-keys.
            return rebuildSubTree1to3(kdNodes, histogram, p);

        } else {
            // The subtree contains more than 3 nodes, so
            // call KdTree.createKdTree to rebuild the subtree,
            // which recycles the nodes and hence invalidates
            // the node argument to this rebuildSubTree function.
            KdTree tree = KdTree.createKdTree(kdNodes, executor, maxSubmitDepth, p);

            // Increment the histogram element.
            if (Constants.ENABLE_HISTOGRAMS) {
                incrementHistogram(histogram, kdNodes.length);
            }

            // Return the root of the rebuilt subtree.
            return tree.root;
        }
    }

    /**
     * <p>
     * The {@code rebuildSubTree1to3} method rebuilds a subtree
     * rooted at a {@code KdNode} that contains 3 nodes or fewer.
     * 
     * @param node - the {@code KdNode} instance
     * @param nodeCount - the number of nodes in the subtree
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    private KdNode rebuildSubTree1to3(final KdNode node,
                                      final int nodeCount,
                                      final long[] histogram,
                                      final int p) {

        // Allocate an array of KdNode references, and walk the subtree
        // to assign to each array element a reference to a KdNode instance.
        final KdNode[] kdNodes = new KdNode[nodeCount];
        getSubTree(node, kdNodes);
        return rebuildSubTree1to3(kdNodes, histogram, p);
    }

    /**
     * <p>
     * The {@code rebuildSubTree1to3} method rebuilds a subtree from
     * an array of {@code KdNode}s that contains 3 nodes or fewer.
     * 
     * @param kdNodes - an array of {@code KdNode} instances
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    private KdNode rebuildSubTree1to3(final KdNode[] kdNodes,
                                      final long[] histogram,
                                      final int p)
    {
        // Increment the histogram element. Note: if '&& !isBalanced(node)' is added
        // to the following 'if' statement, fewer calls to the incrementHistogram
        // method will occur. This observation may explain why no subtree-rebuild
        // operations are reported in histogram bins 0 and 1 for the case where
        // Constants.ENABLE_1TO3==false. But the value of Constants.ENABLE_1TO3
        // affects neither the correctness of the k-d tree, nor the numbers of
        // subtree-rebuild operations reported in histogram bins >= 2.
        if (Constants.ENABLE_HISTOGRAMS) {
            incrementHistogram(histogram, kdNodes.length);
        }

        // Rebuild the subtree and return its root.
        return KdTreeNlogn.buildKdTree1to3(kdNodes, 0, kdNodes.length-1, p);
    }

    /**
     * <p>
     * The {@code rebuildSubTreeSkipRoot1to3} method rebuilds a subtree
     * rooted at a {@code KdNode} that contains 3 nodes or fewer, but
     * omits the root of the subtree from the rebuild.
     * 
     * @param node - the {@code KdNode} instance
     * @param nodeCount - the number of nodes in the subtree
     * @param histogram - a {@code long}[] for counting rebalancing operations
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    private KdNode rebuildSubTreeSkipRoot1to3(final KdNode node,
                                              final int nodeCount,
                                              final long[] histogram,
                                              final int p) {

        // Allocate an array of KdNode references, and walk the subtree
        // to assign to each array element a reference to a KdNode instance.
        final KdNode[] kdNodes = new KdNode[nodeCount];
        getSubTreeSkipRoot(node, kdNodes);
        return rebuildSubTree1to3(kdNodes, histogram, p);
    }

    /**
     * <p>
     * The {@code countNodes} method counts the number of nodes
     * in a subtree rooted at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @return the number of nodes
     * </p>
     */
    private int countNodes(final KdNode node) {

        // Count this node.
        int count = 1;

        // Obtain counts recursively.
        if (node.ltChild != null) {
            count += countNodes(node.ltChild);
        }
        if (node.gtChild != null) {
            count += countNodes(node.gtChild);
        }

        return count;
    }

    /**
     * <p>
     * The {@code getSubTree} method appends to an array the {@code KdNode}
     * references from a subtree rooted at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @param kdNodes - the array of {@code KdNode} instances
     * </p>
     */
    private void getSubTree(final KdNode node,
                            final KdNode[] kdNodes)
    {
        final int[] index = new int[1]; // a kludge to pass an int argument by reference
        index[0]  = 0;
        getSubTree(node, kdNodes, index);
    }

    /**
     * <p>
     * The {@code getSubTree} method appends to an array the {@code KdNode}
     * references from a subtree rooted at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @param kdNodes - the array of {@code KdNode} instances
     * @param index - an int[] array used as a kludge to pass an int by reference
     * </p>
     */
    private void getSubTree(final KdNode node,
                            final KdNode[] kdNodes,
                            final int[] index) {

        // Append the child nodes recursively.
        if (node.ltChild != null) {
            getSubTree(node.ltChild, kdNodes, index);
        }
        if (node.gtChild != null) {
            getSubTree(node.gtChild, kdNodes, index);
        }

        // Reset this node's fields, including its child
        // pointers, and append it to the kdNodes vector.
        node.height = 1;
        node.ltChild = node.gtChild = null;
        kdNodes[index[0]] = node;
        index[0]++;
    }

    /**
     * <p>
     * The {@code countNodesSkipRoot} method counts the number of nodes
     * in a subtree rooted at a (@code KdNode} instance, but omits the
     * root of the subtree from the count.
     * 
     * @param node - the {@code KdNode} instance
     * @return the number of nodes
     * </p>
     */
    private int countNodesSkipRoot(final KdNode node)
    {
        return countNodesSkipRoot(node, 0);
    }

    /**
     * <p>
     * The {@code countNodesSkipRoot} method counts the number of nodes
     * in a subtree rooted at a (@code KdNode} instance, but omits the
     * root of the subtree from the count.
     * 
     * @param node - the {@code KdNode} instance
     * @param depth - the depth in the subtree
     * @return the number of nodes
     * </p>
     */
    private int countNodesSkipRoot(final KdNode node,
                                   final int depth) {

        // Initialize the count.
        int count = 0;

        // Obtain counts recursively.
        if (node.ltChild != null) {
            count += countNodesSkipRoot(node.ltChild, depth+1);
        }
        if (node.gtChild != null) {
            count += countNodesSkipRoot(node.gtChild, depth+1);
        }

        // Count this node unless it's the root.
        if (depth > 0) {
            ++count;
        }

        return count;
    }

    /**
     * <p>
     * The {@code getSubTreeSkipRoot} method appends to an array the {@code KdNode}
     * references from a subtree rooted at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @param kdNodes - the array of {@code KdNode} instances
     * </p>
     */
    private void getSubTreeSkipRoot(final KdNode node,
                                    final KdNode[] kdNodes) {

        final int[] index = new int[1]; // a kludge to pass an int argument by reference
        index[0]  = 0;
        getSubTreeSkipRoot(node, kdNodes, index, 0);
    }

    /**
     * <p>
     * The {@code getSubTreeSkipRoot} method appends to an array the {@code KdNode}
     * references from a subtree rooted at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @param kdNodes - the array of {@code KdNode} instances
     * @param index - an int[] array used as a kludge to pass an int by reference
     * @param depth - the depth in the subtree
     * </p>
     */
    private void getSubTreeSkipRoot(final KdNode node,
                                    final KdNode[] kdNodes,
                                    final int[] index,
                                    final int depth) {

       // Append the child nodes recursively.
        if (node.ltChild != null) {
            getSubTreeSkipRoot(node.ltChild, kdNodes, index, depth+1);
        }
        if (node.gtChild != null) {
            getSubTreeSkipRoot(node.gtChild, kdNodes, index, depth+1);
        }

        // Reset this node's fields, including its child
        // pointers, and append it to the kdNodes vector
        // unless it's the root node.
       if (depth > 0) {
            node.height = 1;
            node.ltChild = node.gtChild = null;
            kdNodes[index[0]] = node;
            index[0]++;
        }
    }

    /**
     * <p>
     * The {@code getSortedTree} method counts and appends the pairs
     * from a {@code KdTree} to a pre-sized vector. Because the tree
     * is traversed in order, the pairs are sorted by their tuples.
     *
     * Calling parameters:
     * 
     * @param node - a KdNode
     * @param coordinates - a Pair<key, value>[]
     * @param index - an int[] that permits an index counter to
     *                be passed by value
     * @return the number of pairs appended
     */
    protected int getSortedPairs(final KdNode node,
                                 final Pair[] coordinates,
                                 int[] index)
    {

        // Initialize the count.
        int count = 0;

        // Obtain counts from < child nodes and append pairs recursively.
        if (node.ltChild != null) {
            count += getSortedPairs(node.ltChild, coordinates, index);
        }

        // Count this node and over-write an entry in the coordinates
        // array by a pair composed of copies of the node's tuple and
        //  the first of its values. Make copies to avoid unintended
        // side effects that could be caused by over-writing randomized
        // coordinate-array elements by sorted-order pairs.
        ++count;
        long[] tuple = new long[node.tuple.length];
        System.arraycopy(node.tuple, 0, tuple, 0, tuple.length);
        String value = new String(node.values.first());
        coordinates[index[0]].putKey(tuple);
        coordinates[index[0]].putValue(value);
        ++index[0];

        // Obtain counts from the > child nodes and append pairs recursively.
        if (node.gtChild != null) {
            count += getSortedPairs(node.gtChild, coordinates, index);
        }

        return count;
    }

    /**
     * <p>
     * The {@code isBalanced} method checks whether the subtree
     * rooted at a (@code KdNode} instance is balanced.
     * 
     * @param node - the {@code KdNode} instance
     * @return true if the subtree is balanced; otherwise, false
     * </p>
     */
    protected static boolean isBalanced(final KdNode node)
    {
        // Get and order the heights at the child nodes.
        final int ltHeight = getHeight(node.ltChild);
        final int gtHeight = getHeight(node.gtChild);
        final int loHeight = (ltHeight < gtHeight) ? ltHeight : gtHeight;
        final int hiHeight = (ltHeight < gtHeight) ? gtHeight : ltHeight;

        if (Constants.AVL_BALANCE)
        {
            // AVL balancing.
            if ( (hiHeight - loHeight) > Constants.HEIGHT_DIFF ) {
                return false;
            } else {
                return true;
            }
        }
        else
        {
            // Red-black balancing. Does the low child exist?
            if (loHeight == 0) {
                // No, the low child does not exist, so test the
                // balance as for a standard AVL tree wherein the
                // difference in child heights must not exceed
                // HEIGHT_DIFF, for which the default value is 1.
                if (hiHeight > Constants.HEIGHT_DIFF) {
                    return false;
                } else {
                    return true;
                }
            } else {
                // Yes, the low child exists, so test the balance
                // as for a red-black tree wherein the height at
                // the high child may be up to twice the height at
                // the low child.
                if ( hiHeight > (loHeight << 1) ) {
                    return false;
                } else {
                    return true;
                }
            }
        }
    }

    /**
     * <p>
     * The {@code computeHeight} method computes the maximum height at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @return the computed height
     * </p>
     */
    protected static int computeHeight(final KdNode node)
    {
        return 1 +  ( (getHeight(node.ltChild) > getHeight(node.gtChild) )
                     ? getHeight(node.ltChild) : getHeight(node.gtChild) ); 
    }

    /**
     * <p>
     * The {@code getHeight} method returns the sub-tree height at a (@code KdNode} instance.
     * 
     * @param node - the {@code KdNode} instance
     * @return the height
     * </p>
     */
    protected static int getHeight(final KdNode node)
    {
        if (node == null) {
            return 0;
        } else {
            return node.height;
        }
    }
    
    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @param verifyLinks - if true, references are checked between AVL nodes and k-d nodes
     * @return the number of nodes in the k-d tree
     */
    protected long verifyKdTree(final boolean verifyLinks)
    {
        // Verify that the doubly linked list has the same number of k-d nodes as the tree.
        if (listSize() != nodeCount) {
            throw new RuntimeException("\n\nlist size = " + listSize() + "  !=  node count =" +
                                       nodeCount + " in KdTreeDynamic.verifyKdTree\n");
        }

        return super.verifyKdTree(verifyLinks);
    }

    /**
     * <p>
     * The {@code verifyKdTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @return the number of nodes in the k-d tree
     */
    protected long verifyKdTree()
    {
        // Verify that the doubly linked list has the same number of k-d nodes as the tree.
        if (listSize() != nodeCount) {
            throw new RuntimeException("\n\nlist size = " + listSize() + "  !=  node count =" +
                                       nodeCount + " in KdTreeDynamic.verifyKdTree\n");
        }

        return super.verifyKdTree();
    }

} // class KdTreeDynamic
