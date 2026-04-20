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

/**
 * <p>
 * The {@code KdTreeLogarithmic} class extends the (@code KdTreeDynamic} class.
 * </p>
 */
public class KdTreeLogarithmic extends KdTreeDynamic {

    private KdTreeDynamic[] kdTrees;
    private AvlTree avlTree;
    protected long[] insertionHistogramLog, deletionHistogramLog;

    KdTreeLogarithmic(final int numDimensions,
                      final ExecutorService executor,
                      final int maxSubmitDepth)
    {
        super(numDimensions, executor, maxSubmitDepth);

        // Create an array to manage the logarithmic k-d trees;
        // by default, each array element is initialized to null.
        kdTrees = new KdTreeDynamic[Constants.MAX_POWER_OF_2];

        // Create histograms to count tree-building operations;
        // by default, each array element is initialized to 0L.
        insertionHistogramDyn = new long[Constants.MAX_POWER_OF_2];
        deletionHistogramDyn = new long[Constants.MAX_POWER_OF_2];

        // Create an instance of AvlTree.
        avlTree = new AvlTree();
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
        if (Constants.ENABLE_DEBUG && coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null in" +
                                       " KdTreeLogarithmic.insert\n");
        }

        // Insert the key-value pair into the AVL Tree.
        boolean inserted = avlTree.insert(coordinate);
        if (Constants.ENABLE_DEBUG && !inserted) {
            throw new RuntimeException("\n\nfailed to insert key-value" +
                                       " pair into AVL tree" +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // If insertion did not create a new AVL node, the key was a duplicate;
        // hence, the value was added to a k-d node's values set in a k-d tree,
        // so no further processing is necessary.
        if (avlTree.insertedNode == null) {
            return true;
        }

        // Insertion created a new AVL node, so the key-value pair will be
        // inserted into a k-d tree in the manner discussed by Overmars and
        // van Leeuwen in "Two general methods for dynamizing decomposable
        // searching problems," Computing 26:155-166, 1981. In particular,
        // the size |B_i| of the ith k-d tree B_i is constrained to lie
        // in the range
        // 
        //     2^(i-2) < |B_i| <= 2^i
        //
        // unless the tree is empty, where the allowed size |B_i| = 0 .
        //
        // The above constraints limit tree sizes as follows (see the
        // analysis further below of insertion into empty B_i).
        //
        // |B_3| = 0 or 2 < |B_3| <= 8
        // |B_2| = 0 or 1 < |B_2| <= 4
        // 0 <= |B_1| <= 2
        // 0 <= |B_0| <= 1
        // 
        // For B_0 and B_1, 2^(i-2) is undefined so a lower limit of 0
        // instead of 2^(i-2) is chosen for these trees.
        //
        // In a departure from the above-referenced discussion by Overmars
        // and van Leeuwen, it is not necessary to build a new k-d tree to
        // insert the key-value pair; instead, the key-value pair is inserted
        // into a non-full (aka sparse) tree that statisfies the constraints
        // 
        //     2^(i-2) < 1 + |B_i| <= 2^i
        //
        // where 1 is the size of the key-value pair and |B_i| is the size
        // of the sparsel tree prior to insertion. Subtracting 1 from all
        // three terms of this equation yields
        //
        //     2^(i-2) <= |B_i| < 2^i
        // 
        // that limits pre-insertion tree sizes as follows.
        //
        // 2 <= |B_3| < 8
        // 1 <= |B_2| < 4
        // 0 <= |B_1| < 2
        // 0 <= |B_0| < 1
        //
        // Inspect the k-d trees in an attempt to find the smallest tree
        // (excluding B_1) that satisfies 2^(i-2) <= |B_i| < 2^i .
        if (Constants.ENABLE_SPARSE_INSERTION) {
            int sparseTree = -1;
            if (getSize(kdTrees[0]) < 1) {
                sparseTree = 0;
            } else if (getSize(kdTrees[1]) < 2) {
                sparseTree = 1;
            } else {
                // Initialize power-of-2 limits for i = 2 and exclude B_0 and B_1
                // from the following 'for' loop.
                int powerOf2Lower = 2, powerOf2Upper = 4;
                for (int i = 2; i < Constants.MAX_POWER_OF_2; ++i) {
                    if (powerOf2Lower <= getSize(kdTrees[i]) &&
                        getSize(kdTrees[i]) < powerOf2Upper) {
                            sparseTree = i;
                            break;
                        }
                    powerOf2Lower <<= 1;
                    powerOf2Upper <<= 1;
                }
            }

            // Was a non-full tree found wherein |2^(i-2) <= |B_i| < 2^i ? 
            if (sparseTree >= 0) {
                // Yes, so insert the key-value pair into that tree, which may be null
                // for B_0 because the getSize method returns 0 for a null or empty tree.
                //
                // Insertion into the non-full tree, and rebalancing that tree,
                // are performed as described by Brown in "A dynamic, self-balancing
                // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
                if (kdTrees[sparseTree] == null) {
                    kdTrees[sparseTree] = new KdTreeDynamic(coordinate.getKey().length,
                                                            executor,
                                                            maxSubmitDepth);
                }
                inserted = kdTrees[sparseTree].insert(coordinate);
                if (Constants.ENABLE_DEBUG && !inserted) {
                    throw new RuntimeException("\n\nfailure to insert coordinate" +
                                                " into sparse tree" + sparseTree +
                                                " in KdTreeLogarithmic.insert\n");
                }
                if (kdTrees[sparseTree].insertedNode == null) {
                    throw new RuntimeException("\n\nnull inserted k-d node" +
                                               " in sparse tree" + sparseTree +
                                               " of KdTreeLogarithmic.insert\n");
                }

                // Create links between the new AVL node and the new k-d node,
                // and return because there is no need to build a new k-d tree.
                avlTree.insertedNode.kdTreeNode = kdTrees[sparseTree].insertedNode;
                kdTrees[sparseTree].insertedNode.avlTreeNode = avlTree.insertedNode;

                // Verify that the size of the tree does not exceed its limits
                // and NOTE that B_0 and B_1 have specific verification tests
                // because 2^(i-2) is undefined for i < 2 so for B_0 and B_1
                // the constraint is 0 <= |B_i| <= 2^i .
                if ( Constants.ENABLE_DEBUG &&
                     getSize(kdTrees[sparseTree]) <= 0 || // |B_i| > 0 after insertion
                     ( (sparseTree == 0 && getSize(kdTrees[0]) > 1) ||
                       (sparseTree == 1 && getSize(kdTrees[1]) > 2) ||
                       (sparseTree >= 2 &&
                         ( (1 << (sparseTree - 2) >= getSize(kdTrees[sparseTree])) ||
                           (getSize(kdTrees[sparseTree]) > (1 << sparseTree)))))) {
                    throw new RuntimeException("\n\ntree " + sparseTree +" size = " +
                                               getSize(kdTrees[sparseTree]) + " in simple" +
                                               " insert of KdTreeLogarithmic.insert\n");
                }

                // Walk the list of k-d tree nodes and ensure that each corresponding
                // AVL tree node indexes the newly created k-d tree in the kdTrees array.
                KdNode node = kdTrees[sparseTree].head;
                while (node != null) {
                    node.avlTreeNode.kdTreeIndex = (short) sparseTree;
                 }
            }
            return true;
        }

        // Insertion of the new key-value pair into tree wherein |2^(i-2) <= |B_i| < 2^i
        // was not requested or else such a tree was not found, so find the smallest tree
        // that is empty.
        int treeIndex; // The index of the replacement tree in kdTrees[treeIndex]
        KdTreeDynamic tree;  // This tree that replaces kdTrees[treeIndex]
        if (getSize(kdTrees[0]) <= 0) {
            // B_0 is empty, so specify B_0 as the smallest tree via treeIndex,
            // and insert the new key-value pair into the smallest tree, B_0.
            treeIndex = 0;
            tree = kdTrees[0];
            if (tree == null) {
                tree = new KdTreeDynamic(coordinate.getKey().length,
                                         executor,
                                         maxSubmitDepth);
            }
            inserted = tree.insert(coordinate);
            if (Constants.ENABLE_DEBUG && !inserted) {
                throw new RuntimeException("\n\nfailure to insert coordinate" +
                                           " into tree 0" +
                                           " in KdTreeLogarithmic.insert\n");
            }
            if (tree.insertedNode == null) {
                throw new RuntimeException("\n\nnull inserted k-d node" +
                                           " in tree 0" +
                                           " of KdTreeLogarithmic.insert\n");
            }
        } else if (getSize(kdTrees[1]) <= 0) {
            // B_0 is not empty but B_1 is empty, set treeIndex to 1,
            // swap B_0 with B1, set kdTreeIndex to reference B_1, and
            // insert the key-value pair into the smallest tree, B_1.
            treeIndex = 1;
            final KdTreeDynamic tmpTree = kdTrees[1];
            kdTrees[1] = kdTrees[0];
            kdTrees[0] = tmpTree;
            kdTrees[1].root.avlTreeNode.kdTreeIndex = 1;
            tree = kdTrees[1]; 
            inserted = tree.insert(coordinate);
             if (Constants.ENABLE_DEBUG && !inserted) {
                throw new RuntimeException("\n\nfailure to insert coordinate" +
                                           " into tree 1" +
                                           " in KdTreeLogarithmic.insert\n");
            }
            if (tree.insertedNode == null) {
                throw new RuntimeException("\n\nnull inserted k-d node" +
                                           " in tree 1" +
                                           " of KdTreeLogarithmic.insert\n");
            }
        } else {
            // Neither B_0 nor B_1 is empty, so scan the k-d trees, beginning
            // with B_2, to find the smallest tree that is empty. Calculate
            // the sum of tree sizes |B_(i-1)| including |B_0| and |B_1|.
            int emptyTree = -1;
            // treeSizeSum includes the size of the to-be-inserted key-value
            // pair (1) as well as |B_0| and |B_1|.
            int treeSizeSum = 1 + getSize(kdTrees[0]) + getSize(kdTrees[1]);
            int powerOf2 = 2; // Initialized to 2^(i-1) for i loop index below.
            for (int i = 2; i < Constants.MAX_POWER_OF_2; ++i) {
                if (getSize(kdTrees[i]) == 0) {
                    emptyTree = i;
                    break;
                }
                treeSizeSum += getSize(kdTrees[i]);
                powerOf2 <<= 1;
            }

            // If no empty tree was found, signal an error.
            if (emptyTree < 0) {
                throw new RuntimeException("\n\nfailed to find empty tree" +
                                           " in KdTreeLogarithmic.insert\n");
            }

            // An empty tree |B_i| was found, so perform the following comparison
            // for i >= 2.
            //
            //     |B_0| + |B_1| + ... + |B_(i-1)| + 1 > 2^(i-1)
            //
            // If the above comparison is true, build a B_i from the key-value
            // pair plus the contents of all trees up to and including B_(i-1).
            // If the above comparison is false, build a B_(i-1) instead.
            // treeSizeSum includes the inserted key-value pair (see above).
            treeIndex = (treeSizeSum > powerOf2) ? emptyTree : emptyTree - 1;

            // Create a new k-d tree with a single k-d node for the inserted
            // key-value pair.
            tree = new KdTreeDynamic(coordinate.getKey().length,
                                     executor,
                                     maxSubmitDepth);

            inserted = tree.insert(coordinate);
            if (Constants.ENABLE_DEBUG && !inserted) {
                throw new RuntimeException("\n\nfailed to insert coordinate into" +
                                           " nascent tree in KdTreeLogarithmic.insert\n");
            }
            if (tree.insertedNode == null) {
                throw new RuntimeException("\n\ninserted node is null" +
                                           " in nascent tree" +
                                           " of KdTreeLogarithmic.insert\n");
            }

            // Append the k-d nodes from trees B_0 ... B_(i-1) to the node list of
            // the newly created k-d tree. Destroy trees B_0 ... B_(i-1) by setting
            // their kdTrees[] elements to null. B_0 and B_1 could be null or empty
            // because they were not included in the search for the first empty tree.
            // But B_2 ... B_(i-1) can be neither null nor empty because the getSize
            // method returned a non-zero result for each of these trees.
            for (int i = 0; i < emptyTree; ++i) {
                if (kdTrees[i] != null) {
                    tree.addList(kdTrees[i]);
                    kdTrees[i] = null;
                }
            }

            // Replace the new k-d tree with a k-d tree created from the list of k-d nodes.
            if (Constants.ENABLE_1TO3 && tree.listNodeCount <= 3) {
                tree = KdTreeDynamic.createKdTree1to3(tree,
                                                      executor,
                                                      maxSubmitDepth,
                                                      insertionHistogramLog);
            } else {
                tree = KdTreeDynamic.createKdTree(tree,
                                                  executor,
                                                  maxSubmitDepth,
                                                  insertionHistogramLog);
            }

            // Verify that the k-d tree contains the correct number of nodes; 
            if (Constants.ENABLE_DEBUG && getSize(tree) != treeSizeSum) {
                throw new RuntimeException("\n\ntree size = " + getSize(tree) +
                                           "  != list size = " + treeSizeSum +
                                           " in KdTreeLogarithmic.insert\n");
            }
        }

        // Verify the correct size of the k-d tree.
        if ( Constants.ENABLE_DEBUG &&
             getSize(tree) <= 0 || // |B_i| > 0 after insertion
             (treeIndex == 0 && getSize(tree) > 1) ||
             (treeIndex == 1 && getSize(tree) > 2) ||
             (treeIndex >= 2 && (1 << (treeIndex - 2) >= getSize(tree) ||
                                 getSize(tree) > (1 << treeIndex)))) {
            throw new RuntimeException("\n\ntree " + treeIndex + " size = " +
                                       getSize(tree) + "in elaborate insert" +
                                       " of KdTreeLogarithmic.insert\n");
        }

        // Create links between the new AVL node and the inserted k-d node.
        avlTree.insertedNode.kdTreeNode = tree.insertedNode;
        tree.insertedNode.avlTreeNode = avlTree.insertedNode;

        // Walk the list of k-d tree nodes and ensure that each corresponding
        // AVL tree node indexes the correct k-d tree in the kdTrees array.
        KdNode node = tree.head;
        while (node != null) {
            node.avlTreeNode.kdTreeIndex = (short) treeIndex;
        }

        // Reference the new k-d tree from the kdTrees array.
        kdTrees[treeIndex] = tree;

        return true;
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
        if (Constants.ENABLE_DEBUG && coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null" +
                                       " in KdTreeLogarithmic.erase\n");
        }

        // Attempt to erase the key-value pair from the AVL tree,
        // and return false if the key was not in the tree.
        if (avlTree.erase(coordinate) == false) {
            return false;
        }

        // If erasure did not remove an AVL node, the key mapped to more than
        // one value; hence, a value was removed from a k-d node's values set
        // but no k-d node was removed from its k-d tree. No further processing
        // is necessary.
        if (avlTree.removedNode == null) {
            return true;
        }

        // Erasure removed an AVL node, so verify the the k-d tree index
        // and erase the corresponding k-d node from its k-d tree.
        //
        // Erasure from the k-d tree, and rebalancing that tree, are performed
        // as described by Brown in "A dynamic, self-balancing k-d tree,"
        // www.arxiv.org/abs/2509.08148, 1-19, 2026.
        final short kdTreeIndex = avlTree.removedNode.kdTreeIndex;
        if ( Constants.ENABLE_DEBUG &&
             (kdTreeIndex < 0) || (kdTreeIndex >= Constants.MAX_POWER_OF_2) ) {
            throw new RuntimeException("\n\nk-d tree index = " + kdTreeIndex +
                                       " is out of range in KdTreeLOgarithmic.erase\n");
        }
        if (kdTrees[kdTreeIndex].erase(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to erase coordinate from" +
                                       " k-d tree in KdTreeLogarithmic.erase\n");
        }
        if (kdTrees[kdTreeIndex].deletedNode == null) {
            throw new RuntimeException("\n\nnull deleted node" +
                                       " in KdTreeLogarithmic.erase\n");
        }

        // Verify that the removed AVL and k-d nodes are correctly linked.
        if ( Constants.ENABLE_DEBUG &&
             ( (avlTree.removedNode.kdTreeNode != kdTrees[kdTreeIndex].deletedNode) ||
               (avlTree.removedNode != kdTrees[kdTreeIndex].deletedNode.avlTreeNode) ) ) {
            throw new RuntimeException("\n\nAVL and k-d nodes are incorrectly" +
                                       " linked in KdTreeLogarithmic.erase\n");
        }

        // The k-d node was removed from its k-d tree, so processing will
        // continue in the manner discussed by Overmars and van Leeuwen
        // in "Two general methods for dynamizing decomposable searching
        // problems," Computing 26:155-166, 1981. In particular, the
        // size |B_i| of the ith k-d tree B_i is constrained to lie in
        // the range
        // 
        //     2^(i-2) < |B_i| <= 2^i
        // 
        // unless the tree is empty, in which case |B_i| = 0.
        // 
        // The above constraints limit tree sizes as follows (see the
        // analysis further below of insertion into empty B_i).
        //
        // |B_3| = 0 or 2 < |B_3| <= 8
        // |B_2| = 0 or 1 < |B_2| <= 4
        // 0 <= |B_1| <= 2
        // 0 <= |B_0| <= 1
        // 
        // For B_0 and B_1, 2^(i-2) is undefined so a lower limit of 0
        // instead of 2^(i-2) is chosen for these trees.
        //
        // A departure from the above-reference discussion by Overmars
        // and van Leeuwen is that deletion of the key-value pair has
        // been followed by rebalancing B_i, as suggested in their article,
        // "Dynamic multi-dimensional data structures based on quad- and
        // k-d trees," Acta Informatica 17: 267-285, 1982.
        //
        // Deletion from the tree, and rebalancing that tree, were
        // performed as described by Brown in "A dynamic, self-balancing
        // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
        //
        // According to Overmars and van Leeuwen, no further processing is
        // necessary if, after deletion of the key-value pair, the following
        // post-deletion constraint is satisfied.
        // 
        //     B_i is empty or 2^(i-2) < |B_i|
        //
        // This constraint omits the constraint |B_i| <= 2^i because
        // the tree-insertion algorithms of the insert method guarantee
        // that this contraint is satisfied.
        //
        // For i < 2, the above post-deletion constraint is analyzed
        // below for the specific cases of B_0 and B_1 because 2^(i-2)
        // is undefined for these cases.
        //
        // Before deletion, B_0 contains 1 k-d node, so deletion of that
        // node produces an empty tree that satisfies the post-deletion
        // constraint. Hence, no further processing is necessary.
        //
        // Before deletion, B_1 contains 2 k-d nodes, so deletion of 1
        // node produces a tree that contains 1 node. This tree violates
        // the post-deletion constraint, so further processing is required.
        // 
        // Consistent with the algorithm below for B_i where i >= 2, and
        // also consistent with the algorithms of the insert method, the
        // following code handles the cases of B_0 and B_1.
        if (kdTreeIndex == 0) {
            // B_0 should be empty due to erasure of its one allowed node,
            // so no further processing is necessary becaise |B_0| = 0
            // is allowed by its constraints. Check that B_0 is indeed empty.
            if (Constants.ENABLE_DEBUG) {
                if (getSize(kdTrees[0]) > 0) {
                    throw new RuntimeException("\n\ntree 0 size = " +
                                               getSize(kdTrees[0]) +
                                               " in empty tree 0 of" +
                                               " KdTreeLogarithmic.erase\n");
                }
            }

            return true;

        } else if (kdTreeIndex == 1 && getSize(kdTrees[1]) < 1) {
            // B_1 is empty, so no further processing is required
            // because its constraints allow |B_1| = 0 .
            //
            // However, becauses B_0 is not empty, it is possible
            // to move B_0's one node to B_1, which is equivalent
            // to the first general case below where |B_(i-1)| > 0
            // and |B_i| = 2^(i-2) but B_1 satifies |B_i| = 0
            // instead of |B_i| = 2^(i-2) .
            //
            // The easiet to move B_0's node to empty B_1 is to
            // swap the trees and then set the kdTreeIndex to
            // reference B_1.
            if (Constants.ENABLE_MOVE_NODE) {
                final KdTreeDynamic tmpTree = kdTrees[1];
                kdTrees[1] = kdTrees[0];
                kdTrees[0] = tmpTree;
                kdTrees[1].root.avlTreeNode.kdTreeIndex = 1;
            }

            return true;

        } else if (getSize(kdTrees[kdTreeIndex]) > (1 << (kdTreeIndex - 2))) {
            // B_i for i >= 2 statifies the contraint |B_1| > 2^(i-2)
            // so no further processing is necessary. But verify the
            // correct size of the k-d tree B_i for i >= 2.
            if (Constants.ENABLE_DEBUG &&
                getSize(kdTrees[kdTreeIndex]) > (1 << kdTreeIndex)) {
                throw new RuntimeException("\n\ntree " + kdTreeIndex +
                                           " size = " +
                                           getSize(kdTrees[kdTreeIndex]) +
                                           " in simple rebuild of" +
                                           " KdTreeLogarithmic.erase\n");
            }

            return true;

        } else {

            // Given the above processing, |B_i| <= 2^(i-2) for i >= 2
            // after deletion, and in fact, for a correctly built tree,
            // after deletion of one node, |B_i| = 2^(i-2) because
            // before deletion |B_i| > 2^(i-2) so declare a new, empty
            // tree that will be subsequently created and then used
            // to build a correctly constrained tree. And declare
            // that tree's index in the kdTrees array.
            KdTreeDynamic tree;
            int treeIndex;

            // Overmars' and van Leeuwen's algorithm determines which
            // k-d tree to build and from which trees to cull the nodes
            // used to build the tree. Destroy the trees from which the
            // nodes are culled by setting their elements to null in the
            // kdTrees array to. The logic below is correct for i >= 2.

            // Is B_(i-1) empty?
            if (getSize(kdTrees[kdTreeIndex-1]) > 0) {
                // B_(i-1) is not empty; hence, the k-d nodes are
                // culled from B_(i-1) and B_i into the new tree.
                tree = new KdTreeDynamic(coordinate.getKey().length,
                                         executor,
                                         maxSubmitDepth);
                tree.addList(kdTrees[kdTreeIndex]);
                tree.addList(kdTrees[kdTreeIndex-1]);
                kdTrees[kdTreeIndex] = kdTrees[kdTreeIndex-1] = null;
                // Is |B_(i-1) > 2^(i-2)?
                if (getSize(kdTrees[kdTreeIndex-1]) > (1 << (kdTreeIndex-2))) {
                    // |B_(i-1)| > 2^(i-2) and |B_i| = 2^(i-2) (see above),
                    // so |B_(i-1)| + |B_i| > 2*2^(i-2) = 2^(i-1) and hence
                    // the combined nodes culled from B_(i-1) and B_i exceed
                    // the post-deletion constraint |B_(i-1)| <= 2^(i-1) so
                    // build B_i using the nodes culled from B_(i-1) and B_i
                    treeIndex = kdTreeIndex;
                } else {
                    // |B_i-1)| + |B_i| <= 2^(i-1) so build B_(i-1)
                    // using the nodes culled from B_(i-1) and B_i
                    treeIndex = kdTreeIndex - 1;
                }
                // Replace either B_i or B_(i-1) with a tree created
                // from the list of nodes.
                if (Constants.ENABLE_1TO3 && tree.listNodeCount <= 3) {
                    kdTrees[treeIndex] =
                        KdTreeDynamic.createKdTree1to3(tree,
                                                       executor,
                                                       maxSubmitDepth,
                                                       insertionHistogramLog);
                } else {
                    kdTrees[treeIndex] =
                        KdTreeDynamic.createKdTree(tree,
                                                   executor,
                                                   maxSubmitDepth,
                                                   insertionHistogramLog);
                }

            // Is B_(i-2) empty?
            } else if (getSize(kdTrees[kdTreeIndex-2]) > 0) {
                // B_(i-1) is empty but B_(i-2) is not empty,
                // so because |B_(i-2)| <= 2^(i-2) and further
                // because |B_i| = 2^(i-2) as discussed above,
                // |B_(i-2)| + |B_i| <= 2*2^(i-2) = 2^(i-1) so
                // build B_(i-1) using the nodes culled from
                // B_i and B_(i-2)
                treeIndex = kdTreeIndex - 1;
                tree = new KdTreeDynamic(coordinate.getKey().length,
                                         executor,
                                         maxSubmitDepth);
                tree.addList(kdTrees[kdTreeIndex]);
                tree.addList(kdTrees[kdTreeIndex-2]);
                kdTrees[kdTreeIndex] = kdTrees[kdTreeIndex-2] = null;
                // Replace B_(i-2) with a tree created from the list of nodes.
                if (Constants.ENABLE_1TO3 && tree.listNodeCount <= 3) {
                    kdTrees[treeIndex] =
                        KdTreeDynamic.createKdTree1to3(tree,
                                                       executor,
                                                       maxSubmitDepth,
                                                       insertionHistogramLog);
                } else {
                    kdTrees[treeIndex] =
                        KdTreeDynamic.createKdTree(tree,
                                                   executor,
                                                   maxSubmitDepth,
                                                   insertionHistogramLog);
                }

            } else {
                // B_(i-1) and B_(i-2) are both empty and, as
                // discussed above, |B_i| = 2^(i-2) so move
                // B_i to B_(i-2) instead of creating a new tree.
                treeIndex = kdTreeIndex - 2;
                tree = kdTrees[treeIndex] = kdTrees[kdTreeIndex];
                kdTrees[kdTreeIndex] = null;
            }

            // Verify the correct size of the k-d tree B_i for i >= 2.
            if ( Constants.ENABLE_DEBUG &&
                getSize(tree) == 0 ||
                (1 << (treeIndex - 2) >= getSize(tree) ||
                getSize(tree) > (1 << treeIndex) ) ) {
                    throw new RuntimeException("\n\ntree " + treeIndex +
                                               " size = " + getSize(tree) +
                                               " in elaborate rebuild of" +
                                               " KdTreeLogarithmic.erase\n");
            }

            // Walk the list of k-d tree nodes and ensure that each corresponding
            // AVL tree node indexes the newly created k-d tree in the kdTrees array.
            KdNode node = tree.head;
            while (node != null) {
                node.avlTreeNode.kdTreeIndex = (short) treeIndex;
            }

            return true;
        }
    }

    /**
     * <p>
     * The {@code verifyTree} method checks that the children of each node of the k-d tree
     * are correctly sorted relative to that node.
     * </p>
     * 
     * @return the number of nodes in the k-d tree
     */
    protected long verifyKdTree()
    {
        // Verify all k-d trees.
        long count = 0L;
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                // Verify that each AVL node contains the correct k-d tree index.
                KdNode node = kdTrees[i].head;
                while (node != null) {
                    if (node.avlTreeNode.kdTreeIndex != (short) i) {
                        throw new RuntimeException("\n\nAVL node index = " + node.avlTreeNode.kdTreeIndex +
                                                "  !=  tree index = " + i +
                                                " in KdTreeLorarithmic.verifyKdTree\n");
                    }
                    node = node.next;
                }
                count += super.verifyKdTree(true);
            }
        }

        return count;
    }

} // class KdTreeLogarithmic
