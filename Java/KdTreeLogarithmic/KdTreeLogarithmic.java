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
            throw new RuntimeException("\n\ncoordinate is null in KdTreeLogarithmic.insert\n");
        }

        // Insert the key-value pair into the AVL Tree.
        boolean inserted = avlTree.insert(coordinate);
        if (Constants.ENABLE_DEBUG && !inserted) {
            throw new RuntimeException("\n\nfailed to insert key-value pair into AVL tree" +
                                       " in KdTreeLogarithmic.insert\n");
        }
        if (avlTree.insertedNode == null) {
            throw new RuntimeException("\n\nnull inserted AVL node" +
                                        " in KdTreeLogarithmic.insert\n");
        }


        // If insertion did not create a new AVL node, the key was a duplicate;
        // hence, the value was added to a k-d node's values set in a k-d tree,
        // so no further processing is necessary.
        if (avlTree.insertedNode != null) {
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
        // unless the tree is empty, in which case the size |B_i| is zero.
        //
        // The above constraints limit tree sizes as follows (see the
        // analysis further below of insertion into empty B_i).
        //
        // B_0 may contain 0 or 1 pair.
        // B_1 may contain 0 or 2 (but not 1) pairs.
        // B_2 may contain 0, 3, or 4 (but not 1 or 2) pairs.
        // B_3 may contain 0, 5, 6, 7, or 8 (but not 1, 2, 3, or 4) pairs.
        //
        // In a departure from the above-referenced discussion by Overmars
        // and van Leeuwen, it is not necessary to build a new k-d tree to
        // insert the key-value pair; instead, the key-value pair is inserted
        // into a non-full tree that statisfies 2^(i-2) < 1 + |B_i| <= 2^i
        // where 1 is the size of the key-value pair and |B_i| is the size
        // of the non-full tree. Subtracting 1 from all three terms of this
        // comparison yields 2^(i-2) <= |B_i| < 2^i so for example, before
        // insertion of a new key-value pair, B_2 may contain 2 or 3 pairs,
        // and B_3 may contain 4, 5, 6, or 7 pairs.
        //
        // Insertion into the non-full tree, and rebalancing that tree,
        // are performed as described by Brown in "A dynamic, self-balancing
        // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
        //
        // Inspect the k-d trees in an attempt to find the smallest tree wherein
        // 2^(i-2) <= |B_i| < 2^i but NOTE that getSize() returns zero for
        // either a null tree or an empty tree. This getSize() behavior creates
        // a special case for B_1 because B_1 may contain 0 or 2 pairs but NOT
        // 1 pair, and a search based on the getSize method for a tree wherein
        // |B_i| < 2^i would identify empty B_1 as such a tree. However, insertion
        // of a new pair nto empty B_1 would cause it to contain 1 pair, which
        // violates the allowed tree size for B_1 discussed above. Moreover, the
        // comparison 2^(i-2) <= |B_i| is not possible for i < 2 because 2^(i-2)
        // cannot be represented as an integer. For these reasons, B_1 is excluded
        // from the 'for' loop below for the tree wherein 2^(i-2) <= |B_i| < 2^i
        // and B_0 is inspected outside of that 'for' loop.
        if (Constants.ENABLE_NON_FULL_INSERTION) {
            int sparseTree = -1;
            if (getSize(kdTrees[0]) == 0)  {
                sparseTree = 0;
            } else {
                 // Initialize power-of-2 limits for i = 2 per above discussion.
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
            // Was a non-full tree found?
            if (sparseTree >= 0) {
                // Yes, so insert the key-value pair into that tree, which may be null.
                if (kdTrees[sparseTree] == null) {
                    kdTrees[sparseTree] = new KdTreeDynamic(coordinate.getKey().length,
                                                            executor,
                                                            maxSubmitDepth);
                }
                inserted = kdTrees[sparseTree].insert(coordinate);
                if (Constants.ENABLE_DEBUG && !inserted) {
                    throw new RuntimeException("\n\nfailure to insert coordinate" +
                                                " in KdTreeLogarithmic.insert\n");
                }
                if (kdTrees[sparseTree].insertedNode == null) {
                    throw new RuntimeException("\n\nnull inserted k-d node" +
                                               " in KdTreeLogarithmic.insert\n");
                }

                // Create links between the new AVL node and the new k-d node,
                // and return because there is no need to build a new k-d tree.
                avlTree.insertedNode.kdTreeNode = kdTrees[sparseTree].insertedNode;
                kdTrees[sparseTree].insertedNode.avlTreeNode = avlTree.insertedNode;

                // Verify that the size of the tree does not exceed its limits
                // and NOTE that B_0 and B_1 have specific verification tests.
                if ( Constants.ENABLE_DEBUG &&
                     getSize(kdTrees[sparseTree]) == 0 ||
                     ( (sparseTree == 0 && getSize(kdTrees[0]) != 1) ||
                       (sparseTree == 1 && getSize(kdTrees[1]) != 2) ||
                       (sparseTree > 1 &&
                         ( (1 << (sparseTree - 2) >= getSize(kdTrees[sparseTree]) ) ||
                           (getSize(kdTrees[sparseTree]) > (1 << sparseTree) ) ) ) ) ) {
                    throw new RuntimeException("\n\ntree " + sparseTree +" size = " +
                                               getSize(kdTrees[sparseTree]) + " in simple" +
                                               " insert of KdTreeLogarithmic.insert\n");
                }

                // Walk the list of k-d tree nodes and ensure that each corresponding
                // AVL tree node indexes the newly created k-d tree in the kdTrees array.
                KdNode node = kdTrees[sparseTree].getHead();
                while (node != null) {
                    node.avlTreeNode.kdTreeIndex = (short) sparseTree;
                 }
            }
            return true;
        }

        // Either insertion of the new key-value pair into a non-empty tree
        // was not requested or a non-empty tree was not found.
        //
        // Scan the k-d trees, beginning with B_1, to find the smallest tree
        // that is empty. Calculate the sum of tree sizes |B_(i-1)| including
        // |B_0|. Do not allow B_0 to be identified as the empty B_i because
        // B_(i-1) might be used instead of B_i below and because B_(-1) does
        // not exist.
        int emptyTree = -1;
        int powerOf2 = 1; // Initialized to 2^(i-1) for i loop index below.
        int treeSizeSum = 1 + getSize(kdTrees[0]); // Include the key-value pair.
        for (int i = 1; i < Constants.MAX_POWER_OF_2; ++i) {
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

        // The empty tree |B_i| was found, so perform the following comparison
        // that is valid for i < 0, and hence the reason that the above search
        // for the smallest empty tree ignores B_0.
        //
        //     |B_0| + |B_1| + ... + |B_(i-1)| + 1 > 2^(i-1)
        //
        // If the above comparison is true, build a B_i from the key-value
        // pair plus the contents of all trees up to and including B_(i-1).
        // If the above comparison is false, build a B_(i-1) instead.
        // treeSizeSum includes the inserted key-value pair (see above).
        //
        // Analysis of the above comparison for empty tree B_1 is useful.
        // If B_0 is empty, the result of the comparison 1 > 2^0 is false,
        // so the new pair will be inserted into B_0. If B_0 is not empty,
        // |B_0| = 1 and the result of the comparison 2 > 2^0 is true, so
        // the new pair will be inserted into B_1 such that |B_1| = 2.
        // This analysis reveals that insertion cannot produce |B_1| = 1,
        // consistent the earlier discussion of the tree sizes permitted
        // by the constraints 2^(i-2) < |B_i| <= 2^i
        int treeIndex = (treeSizeSum > powerOf2) ? emptyTree : emptyTree - 1;

        // Create a new k-d tree with a single k-d node for the inserted key-value pair.
        KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                               executor,
                                               maxSubmitDepth);

        inserted = tree.insert(coordinate);
        if (Constants.ENABLE_DEBUG && !inserted) {
            throw new RuntimeException("\n\nfailed to insert coordinate into" +
                                       " k-d tree in KdTreeLogarithmic.insert\n");
        }
        if (tree.insertedNode == null) {
            throw new RuntimeException("\n\ninserted node is null" +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // Create links between the new AVL node and the single k-d node.
        avlTree.insertedNode.kdTreeNode = tree.insertedNode;
        tree.insertedNode.avlTreeNode = avlTree.insertedNode;

        // Append the k-d nodes from trees B_0 ... B_(i-1) to the node list of
        // the newly created k-d tree. Destroy trees B_0 ... B_(i-1) by setting
        // their kdTrees[] elements to null. B_0 may already be null because it
        // was not included in the search for the first empty tree. However,
        // B_1 ... B_(i-1) can be neither null nor empty because the getSize
        // method returned a non-zero result for each of these trees.
        if (kdTrees[0] != null) {
            tree.addList(kdTrees[0]);
            kdTrees[0] = null;
        }
        for (int i = 1; i < emptyTree; ++i) {
            tree.addList(kdTrees[i]);
            kdTrees[i] = null;
        }

        // Replace the new k-d tree with a k-d tree created from the list of k-d nodes.
        if (Constants.ENABLE_1TO3 && tree.listSize() <= 3) {
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
        // treeSizeSum includes the inserted key-value pair (see above).
        if (Constants.ENABLE_DEBUG && getSize(tree) != treeSizeSum) {
            throw new RuntimeException("\n\ntree size = " + getSize(tree) +
                                       "  != list size = " + treeSizeSum +
                                       " in KdTreeLogarithmic.insert\n");
        }
        if ( Constants.ENABLE_DEBUG &&
             getSize(tree) == 0 ||
             (emptyTree == 0 && getSize(tree) != 1) ||
             (emptyTree == 1 && getSize(tree) != 2) ||
             (emptyTree >= 2 && (1 << (emptyTree - 2) >= getSize(tree) ||
                                 getSize(tree) > (1 << emptyTree) ) ) ) {
            throw new RuntimeException("\n\ntree " + emptyTree + " size = " +
            getSize(tree) + "in elaborate insert of KdTreeLogarithmic.insert\n");
        }

        // Walk the list of k-d tree nodes and ensure that each  corresponding
        // AVL tree node indexes the newly created k-d tree in the kdTrees array.
        KdNode node = tree.getHead();
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
            throw new RuntimeException("\n\ncoordinate is null in KdTreeLogarithmic.erase\n");
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
        // unless the tree is empty, in which case the size is zero.
        // 
        // The above constraints limit tree sizes as follows (see the
        // analysis further below of insertion into empty B_i).
        //
        // B_0 may contain 0 or 1 pair.
        // B_1 may contain 0 or 2 (but not 1) pairs.
        // B_2 may contain 0, 3, or 4 (but not 1 or 2) pairs.
        // B_3 may contain 0, 5, 6, 7, or 8 (but not 1, 2, 3, or 4) pairs.
        //
        // A departure from the above-reference discussion by Overmars
        // and van Leeuwen is that deletion of the key-value pair has
        // been followed by rebalancing B_i, as suggested in their article,
        // "Dynamic multi-dimensional data structures based on quad- and
        // k-d trees," Acta Informatica 17: 267-285, 1982.
        //
        // Deletion from the tree, and rebalancing that tree, have been
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
            // B_0 is empty, so no further processing is necessary.
            return true;
        } else if (kdTreeIndex == 1) {
            // B_1 contains 1 k-d node, so further processing is required
            // to build a correct tree.
           if (getSize(kdTrees[0]) > 0) {
                // B_0 and B_1 each contain 1 node, so B_1's node is
                // added to B_0's node list, the list is used to build
                // a new B_1, and then B_0 is destroyed.
                kdTrees[0].addList(kdTrees[1]);
                KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                                       executor,
                                                       maxSubmitDepth);
                kdTrees[1] = createKdTree(tree,
                                          executor,
                                          kdTreeIndex,
                                          deletionHistogramDyn);
                kdTrees[0] = null;

                // Increment the histogram element.
                if (Constants.ENABLE_HISTOGRAMS) {
                    incrementHistogram(deletionHistogramLog, 2);
                }
            } else {
                // B_0 is empty and B_1 contains 1 node, so B_1 replaces
                // B_0 and then B_1 is destroyed.
                kdTrees[0] = kdTrees[1];
                kdTrees[1] = null;
 
                // Increment the histogram element.
                if (Constants.ENABLE_HISTOGRAMS) {
                    incrementHistogram(deletionHistogramLog, 1);
                }
           }
        } else if (getSize(kdTrees[kdTreeIndex]) > (1 << (kdTreeIndex - 2))) {
            return true;
        }

        // Given the above processing, |B_i| <= 2^(i-2) for i >= 2
        // after deletion, and in fact, for a correctly built tree,
        // |B_i| = 2^(i-2) because before deletion |B_i| > 2^(i-2)
        // so create a new, empty k-d tree that will be used to
        // build a correctly constrained tree. 
        KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                               executor,
                                               maxSubmitDepth);

        // Overmars' and van Leeuwen's algorithm determines which
        // k-d tree to build and from which trees to cull the nodes
        // used to build the tree. Destroy the trees from which the
        // nodes are culled by setting their elements in the kdTrees
        // array to null. The logic below is correct for i >= 2.
        int treeIndex;
        // Is B_(i-1) empty?
        if (getSize(kdTrees[kdTreeIndex-1]) > 0) {
            // B_(i-1) is not empty, so the k-d nodes are
            // culled from B_(i-1) and B_i
            tree.addList(kdTrees[kdTreeIndex-1]);
            tree.addList(kdTrees[kdTreeIndex]);
            kdTrees[kdTreeIndex-1] = kdTrees[kdTreeIndex] = null;
            // Is |B_(i-1) > 2^(i-2)?
            if (getSize(kdTrees[kdTreeIndex-1]) > (1 << (kdTreeIndex-2))) {
                // |B_(i-1)| > 2^(i-2) and |B_i| = 2^(i-2) (see above),
                // so |B_(i-1)| + |B_i| > 2*2^(i-2) = 2^(i-1) hence
                // the combined nodes culled from B_(i-1) and B_i exceed
                // the post-deletion constraint |B_(i-1)| <= 2^(i-1) so
                // build B_i using the nodes culled from B_(i-1) and B_i
                treeIndex = kdTreeIndex;
            } else {
                // |B_i-1)| + |B_i| <= 2^(i-1) so build B_(i-1)
                // using the nodes culled from B_(i-1) and B_i
                treeIndex = kdTreeIndex - 1;
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
            tree.addList(kdTrees[kdTreeIndex]);
            tree.addList(kdTrees[kdTreeIndex-2]);
            kdTrees[kdTreeIndex] = kdTrees[kdTreeIndex-2] = null;
         } else {
            // B_(i-1) and B_(i-2) are both empty and, as
            // discussed above, |B_i| = 2^(k-2) so build
            // B_(i-2) using the nodes culled from B_i
            treeIndex = kdTreeIndex - 2;
            tree.addList(kdTrees[kdTreeIndex]);
            kdTrees[kdTreeIndex] = null;
        }

        // Replace the new k-d tree with a tree created from the list of nodes.
        if (Constants.ENABLE_1TO3 && tree.listSize() <= 3) {
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

        if ( Constants.ENABLE_DEBUG &&
             getSize(tree) == 0 ||
             (treeIndex == 0 && getSize(tree) != 1) ||
             (treeIndex == 1 && getSize(tree) != 2) ||
             (treeIndex >= 2 && (1 << (treeIndex - 2) >= getSize(tree) ||
                                 getSize(tree) > (1 << treeIndex) ) ) ) {
            throw new RuntimeException("\n\ntree " + treeIndex +
                                       " size = " + getSize(tree) +
                                       " in KdTreeLogarithmic.erase\n");
        }

        // Walk the list of k-d tree nodes and ensure that each corresponding
        // AVL tree node indexes the newly created k-d tree in the kdTrees array.
        KdNode node = tree.getHead();
        while (node != null) {
            node.avlTreeNode.kdTreeIndex = (short) treeIndex;
        }

        // Reference the new k-d tree from the kdTrees array.
        kdTrees[treeIndex] = tree;

        return true;
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
                KdNode node = kdTrees[i].getHead();
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
