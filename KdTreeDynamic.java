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

/**
 * <p>
 * @author Russell A. Brown
 */

import java.util.concurrent.ExecutorService;
import java.util.TreeSet;

/**
 * <p>
 * The {@code KdTreeDynamic}  class extends the (@code KdTree} class and
 * accessess KdTree.root so that an instance of KdTreeDynamic
 * may use KdTree member functions that require access to root.
  * </p>
 */
public class KdTreeDynamic extends KdTree {

    private int maxSubmitDepth = -1;
    private int cutoff = Constants.MULTI_THREAD_CUTOFF;
    private boolean inserted = false, erased = false, changed = false;
    private ExecutorService executor;

    // This constructor does not set KdTree.root
    KdTreeDynamic(final int maxSubmitDepth,
                  final int cutoff,
                  final ExecutorService executor)
    {
        this.maxSubmitDepth = maxSubmitDepth;
        this.cutoff = cutoff;
        this.executor = executor;
    }

    // This constructor set KdTree.root
    KdTreeDynamic(final int maxSubmitDepth,
                  final int cutoff,
                  final ExecutorService executor,
                  KdNode node)
    {
        this.maxSubmitDepth = maxSubmitDepth;
        this.cutoff = cutoff;
        this.executor = executor;
        root = node; // the root of the KdTree instance
    }

    protected boolean isEmpty()
    {
        return (root == null);
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
        inserted = changed = false;
        if (root != null) {
           root = insert(root, coordinate.getKey(), coordinate.getValue(), 0);

            // Has the height changed due to an insertion?
            if (changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(root) ) {
                    root.height = computeHeight(root);
                } else {
                    root = rebuildSubTree(root, 0);
                }
            }
        } else {
            root = new KdNode(coordinate.getKey().length);
            root.height = 1;
            System.arraycopy(coordinate.getKey(), 0, root.tuple, 0, root.tuple.length);
            root.values.add(coordinate.getValue());
            inserted = changed = true;
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
     * @param q - the leading dimension that is permuted cyclically
     * @return the (@code KdNode) root of the subtree
     * </p>
     */
    private KdNode insert(final KdNode node,
                         final long[] key,
                         final String value,
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
                nodePtr.ltChild = insert(nodePtr.ltChild, key, value, p+1);
            } else {
                // The tree does not contain the key, so insert it.
                nodePtr.ltChild = new KdNode(dim);
                nodePtr.ltChild.height = 1;
                System.arraycopy(key, 0, nodePtr.ltChild.tuple, 0, dim);
                nodePtr.ltChild.values.add(value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) > 0L) {
            if (nodePtr.gtChild != null) {
                nodePtr.gtChild = insert(nodePtr.gtChild, key, value, p+1);
            } else {
                // The tree does not contain the key, so insert it.
                nodePtr.gtChild = new KdNode(dim);
                nodePtr.gtChild.height = 1;
                System.arraycopy(key, 0, nodePtr.gtChild.tuple, 0, dim);
                nodePtr.gtChild.values.add(value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else {
            // The tree contains the key, so insert the value,
            // although the values set might already contain
            // it. The tree height doesn't change.
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
            if ( isBalanced(nodePtr) ) {
                nodePtr.height = computeHeight(nodePtr);
            } else {
                nodePtr = rebuildSubTree(nodePtr, p);
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
    protected boolean erase(final Pair coordinate) {

        erased = changed = false;
        if (root != null) {
            root = erase(root, coordinate.getKey(), coordinate.getValue(), 0, false);

            // Has the height changed due to an erasure?
            if (root != null && changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(root) ) {
                    root.height = computeHeight(root);
                } else {
                    root = rebuildSubTree(root, 0);
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
     * @param q - the leading dimension that is permuted cyclically
     * @param clearSet - clear all values in the values set
     * @return the (@code KdNode) root of the subtree
     * </p>
     */
    private KdNode erase(final KdNode node,
                         final long[] key,
                         final String value,
                         final int q,
                         final boolean clearSet) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        final int dim = node.tuple.length;
        final int p = (q < dim) ? q : 0;

        // The return value because the node argument to this function is read only.
        KdNode nodePtr = node;

        // Determine which child to search recursively for a deletion point.
        if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) < 0L) {
            if (nodePtr.ltChild != null) {
                nodePtr.ltChild = erase(nodePtr.ltChild, key, value, p+1, clearSet);

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
                    if ( isBalanced(nodePtr) ) {
                        nodePtr.height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, p);
                    }
                }
            } else {
                // The tree does not contain the key.
                erased = changed = false;
            }
        } else if (MergeSort.superKeyCompare(key, nodePtr.tuple, p) > 0L) {
            if (nodePtr.gtChild != null) {
                nodePtr.gtChild = erase(nodePtr.gtChild, key, value, p+1, clearSet);
             
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
                    if ( isBalanced(nodePtr) ) {
                        nodePtr.height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, p);
                    }
                }
            } else {
                // The tree does not contain the key.
                erased = changed = false;
            }
        } else {
            // The tree contains the key, so if clearSet is false,
            // and if the values set does not contain the value, no
            // value is removed and the tree height doesn't change.
            if ( !clearSet && !nodePtr.values.contains(value) ) {
                erased = changed = false;
            } else {
                // Either clearSet is true, or the values set
                // contains the value. So if clearSet is true,
                // clear the values set; otherwise, remove the
                // value from the values set.
                if (clearSet) {
                    nodePtr.values.clear();
                } else {
                    nodePtr.values.remove(value);
                }

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
                        // Yes, the node has only a < child, so does the subtree
                        // rooted at the < child contain <= 3 nodes?
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
                            final KdNode tempPtr = nodePtr;
                            nodePtr = rebuildSubTree1to3(nodePtr.ltChild, nodeCount, p);
                            tempPtr.ltChild = null; // Prevent recursive deletion.
                        }
                        else
                        {
                            // No, the < child subtree contains > 3 nodes. So, find
                            // the immediate predecessor node, copy the tuple and the
                            // values set from that predecessor node to the one-child
                            // node, delete the predecessor node recursively (clearing
                            // its values set), and recompute the height along the path
                            // back to the < child, including that child.
                            KdNode predecessor = nodePtr.ltChild;
                            predecessor = findPredecessor(nodePtr.ltChild, predecessor, p, p+1);
                            System.arraycopy(predecessor.tuple, 0, nodePtr.tuple, 0, dim);
                            nodePtr.values.clear();
                            nodePtr.values.addAll(predecessor.values);
                            // value is a dummy argument because the clearSet argument is true
                            nodePtr.ltChild = erase(nodePtr.ltChild, nodePtr.tuple, value, p+1, true);

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
                            if ( isBalanced(nodePtr) ) {
                                nodePtr.height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, p);
                            }
                        }
                    }
                    // Does the node have only a > child?
                    else if (nodePtr.gtChild != null && nodePtr.ltChild == null) {
                        // Yes, the node has only a > child, so does the subtree
                        // rooted at the > child contain <= 3 nodes?
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
                            final KdNode tempPtr = nodePtr;
                            nodePtr = rebuildSubTree1to3(nodePtr.gtChild, nodeCount, p);
                            tempPtr.gtChild = null; // Prevent recursive deletion.
                        }
                        else
                        {
                            // No, the > child subtree contains > 3 nodes. So, find
                            // the immediate successor node, copy the tuple and the
                            // values set from that successor node to the one-child
                            // node, delete the successor node recursively (clearing
                            // its values set), and recompute the height along the path
                            // back to the > child, including that child.
                            KdNode successor = nodePtr.gtChild;
                            successor = findSuccessor(nodePtr.gtChild, successor, p, p+1);
                            System.arraycopy(successor.tuple, 0, nodePtr.tuple, 0, dim);
                            nodePtr.values.clear();
                            nodePtr.values.addAll(successor.values);
                            // value is a dummy argument because the clearSet argument is true
                            nodePtr.gtChild = erase(nodePtr.gtChild, nodePtr.tuple, value, p+1, true);

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
                            if ( isBalanced(nodePtr) ) {
                                nodePtr.height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, p);
                            }
                        }
                    }
                    // If the node has no children, delete the node.
                    else if (nodePtr.ltChild == null && nodePtr.gtChild == null) {
                        nodePtr = null; // Garbage collection will occur anyway without this.
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
                            KdNode tempPtr = nodePtr;
                            nodePtr = rebuildSubTreeSkipRoot1to3(nodePtr, nodeCount, p);
                            tempPtr.ltChild = tempPtr.gtChild = null; // Prevent recursive deletion.
                        }
                        else
                        {
                            // No, the subtree rooted at this two-child node contains
                            // > 3 nodes, excluding this two-child node, so replace this
                            // two-child node by either its predecessor or its successor.
                            //
                            // If Constants.ENABLE_PREFERRED_TEST is true, select the replacement
                            // node from the taller of the child subtrees.
                            if (Constants.ENABLE_PREFERRED_NODE
                                && (getHeight(nodePtr.ltChild) >= getHeight(nodePtr.gtChild)))
                            {
                                // Find the node with the largest super-key in the
                                // subtree rooted at the < child, which is the
                                // predecessor node. Copy the predecessor node's tuple
                                // and values set to this two-child node, delete the
                                // predecessor node recursively (clearing its values set),
                                // and recompute the heights along the path from the
                                // predecessor node to (but excluding) this two-child node.
                                KdNode predecessor = nodePtr.ltChild;
                                predecessor = findPredecessor(nodePtr.ltChild, predecessor, p, p+1);
                                System.arraycopy(predecessor.tuple, 0, nodePtr.tuple, 0, dim);
                                nodePtr.values.clear();
                                nodePtr.values.addAll(predecessor.values);
                                // value is a dummy argument because the clearSet argument is true
                                nodePtr.ltChild = erase(nodePtr.ltChild, nodePtr.tuple, value, p+1, true);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( isBalanced(nodePtr) ) {
                                    nodePtr.height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, p);
                                }
                            }
                            else
                            {
                                // Find the node with the smallest super-key in the
                                // subtree rooted at the > child, which is the
                                // successor node. Copy the successor node's tuple
                                // and values set to this two-child node, delete the
                                // successor node recursively (clearing its values set),
                                // and recompute the heights along the path from the
                                // successor node to (but excluding) this two-child node.
                                KdNode successor = nodePtr.gtChild;
                                successor = findSuccessor(nodePtr.gtChild, successor, p, p+1);
                                System.arraycopy(successor.tuple, 0, nodePtr.tuple, 0, dim);
                                nodePtr.values.clear();
                                nodePtr.values.addAll(successor.values);
                                // value is a dummy argument because the clearSet argument is true
                                nodePtr.gtChild = erase(nodePtr.gtChild, nodePtr.tuple, value, p+1, true);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( isBalanced(nodePtr) ) {
                                    nodePtr.height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, p);
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
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    protected KdNode rebuildSubTree(final KdNode node,
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
            return rebuildSubTree1to3(node, kdNodes, p);

        } else {
            // The subtree contains more than 3 nodes, so
            // call KdTree.createKdTree to rebuild the subtree,
            // which recycles the nodes and hence invalidates
            // the node argument to this rebuildSubTree function.
            //
            // Do not use mutiple threads to build the subtree, even
            // though multiple threads are available, unless the size of the
            // subtree is sufficiently large to justify spawning child threads.
            KdTree tree;
            if (count < cutoff) {
                tree = KdTree.createKdTree(kdNodes, executor, -1, p);
            } else {
                tree = KdTree.createKdTree(kdNodes, executor, maxSubmitDepth, p);
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
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    private static KdNode rebuildSubTree1to3(final KdNode node,
                                             final int nodeCount,
                                             final int p) {

        // Allocate an array of KdNode references, and walk the subtree
        // to assign to each array element a reference to a KdNode instance.
        final KdNode[] kdNodes = new KdNode[nodeCount];
        getSubTree(node, kdNodes);
        return rebuildSubTree1to3(node, kdNodes, p);
    }

    private static KdNode rebuildSubTree1to3(final KdNode node,
                                             final KdNode[] kdNodes,
                                             final int p)
    {

        // The return value because the node argument to this function is read only.
        KdNode ptr = node;

        if (kdNodes.length == 1) {

            // The subtree contains 1 node, so it is the root of the rebuilt subtree.
            ptr = kdNodes[0];
            ptr.height = 1;

        } else if (kdNodes.length == 2) {

            // The subtree contains 2 nodes, so the first node is the root of the rebuilt
            // subtree. Compare super-keys to determine whether second node is the < child
            // or the > child of the first node.
            ptr = kdNodes[0];
            ptr.height = 2;
            if (MergeSort.superKeyCompare(kdNodes[0].tuple, kdNodes[1].tuple, p) < 0L) {
                ptr.gtChild = kdNodes[1];
                ptr.gtChild.height = 1;
            } else {
                ptr.ltChild = kdNodes[1];
                ptr.ltChild.height = 1;
            }

        } else if (kdNodes.length == 3) {

            // The subtree contains 3 nodes, so compare their super-keys to determine which
            // node is the median node. The node with the smallest super-key is the < child.
            // The node with the larges super-key is the > child.
            if (MergeSort.superKeyCompare(kdNodes[0].tuple, kdNodes[1].tuple, p) < 0L) {
                // kdNodes[0].tuple < kdNodes[1].tuple
                if (MergeSort.superKeyCompare(kdNodes[1].tuple, kdNodes[2].tuple, p) < 0L) {
                    // kdNodes[0].tuple < kdNodes[1].tuple < kdNodes[2].tuple
                    ptr = kdNodes[1];
                    ptr.ltChild = kdNodes[0];
                    ptr.gtChild = kdNodes[2];
                } else {
                    // kdNodes[0].tuple < kdNodes[1].tuple; kdNodes[2].tuple < kdNodes[1].tuple
                    if (MergeSort.superKeyCompare(kdNodes[0].tuple, kdNodes[2].tuple, p) < 0L) {
                        // kdNodes[0].tuple < kdNodes[2].tuple < kdNodes[1].tuple
                        ptr = kdNodes[2];
                        ptr.ltChild = kdNodes[0];
                        ptr.gtChild = kdNodes[1];
                    } else {
                        // kdNodes[2].tuple < kdNodes[0].tuple < kdNodes[1].tuple
                        ptr = kdNodes[0];
                        ptr.ltChild = kdNodes[2];
                        ptr.gtChild = kdNodes[1];
                    }
                }
            } else {
                // kdNodes[1].tuple < kdNodes[0].tuple
                if (MergeSort.superKeyCompare(kdNodes[0].tuple, kdNodes[2].tuple, p) < 0L) {
                    // kdNodes[1].tuple < kdNodes[0].tuple < kdNodes[2].tuple
                    ptr = kdNodes[0];
                    ptr.ltChild = kdNodes[1];
                    ptr.gtChild = kdNodes[2];
                } else {
                    // kdNodes[1].tuple < kdNodes[0].tuple; kdNodes[2].tuple < kdNodes[0].tuple
                    if (MergeSort.superKeyCompare(kdNodes[1].tuple, kdNodes[2].tuple, p) < 0L) {
                        // kdNodes[1].tuple < kdNodes[2].tuple < kdNodes[0].tuple
                        ptr = kdNodes[2];
                        ptr.ltChild = kdNodes[1];
                        ptr.gtChild = kdNodes[0];
                    }
                    else {
                        // kdNodes[2].tuple < kdNodes[1].tuple < kdNodes[0].tuple
                        ptr = kdNodes[1];
                        ptr.ltChild = kdNodes[2];
                        ptr.gtChild = kdNodes[0];
                    }
                }
            }
            ptr.ltChild.height = ptr.gtChild.height = 1;
            ptr.height = 2;

        } else {

            // This is an illegal condition that should never occur.
            throw new RuntimeException("\n\n" + kdNodes.length + " KdNode instances in rebuildSubTree1to3\n");
        }

        return ptr;
    }

    /**
     * <p>
     * The {@code rebuildSubTreeSkipRoot1to3} method rebuilds a subtree
     * rooted at a {@code KdNode} that contains 3 nodes or fewer, but
     * omits the root of the subtree from the rebuild.
     * 
     * @param node - the {@code KdNode} instance
     * @param nodeCount - the number of nodes in the subtree
     * @param p - the leading dimension
     * @return the root {@code KdNode} of the rebalanced subtree
     * </p>
     */
    private static KdNode rebuildSubTreeSkipRoot1to3(final KdNode node,
                                                     final int nodeCount,
                                                     final int p) {

        // Allocate an array of KdNode references, and walk the subtree
        // to assign to each array element a reference to a KdNode instance.
        final KdNode[] kdNodes = new KdNode[nodeCount];
        getSubTreeSkipRoot(node, kdNodes);
        return rebuildSubTree1to3(node, kdNodes, p);
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
    private static int countNodes(final KdNode node) {

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
    private static void getSubTree(final KdNode node,
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
    private static void getSubTree(final KdNode node,
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
    private static int countNodesSkipRoot(final KdNode node)
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
    private static int countNodesSkipRoot(final KdNode node,
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
    private static void getSubTreeSkipRoot(final KdNode node,
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
    private static void getSubTreeSkipRoot(final KdNode node,
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
     * The {@code getHeight} method returns the subtree height at a (@code KdNode} instance.
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
} // class KdTreeDynamic
