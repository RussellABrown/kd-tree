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
        if (Constants.ENABLE_DEBUG && avlTree.insert(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to insert key-value pair into AVL tree" +
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
        // unless the tree is empty, in which case the size is zero.
        //
        // In a departure from the above-referenced discussion by Overmars
        // and van Leeuwen, a new k-d tree is not perforce built to insert
        // the key-value pair; instead, the key-value pair is inserted into
        // a non-empty tree if the tree statisfies the condition |B_i| < 2^i
        // that is the same condition as |B_i| + 1 <= 2^i where 1 represents
        // the size of the new k-d node that holds the key-value pair.
        //
        // Inspect the k-d trees, beginning with the smallest tree, to find
        // the smallest tree for which |B_i| < 2^i while using the treeSize
        // method to obtain the tree's size, because treeSize returns zero
        // for either a null tree or an empty tree.
        if (Constants.ENABLE_NON_FULL_INSERTION) {
            int firstNonFullTree = -1;
            int powerOf2 = 1; // Initialized to 2^i for i computed below
            for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
                if (treeSize(kdTrees[i]) < powerOf2) {
                    firstNonFullTree = i;
                    break;
                }
                powerOf2 <<= 1;
            }
            // Was a non-full tree found?
            if (firstNonFullTree >= 0) {
                // Yes, so insert the key-value pair into that tree, which may be null.
                if (kdTrees[firstNonFullTree] == null) {
                    kdTrees[firstNonFullTree] = new KdTreeDynamic(coordinate.getKey().length,
                                                                  executor,
                                                                  maxSubmitDepth);
                }
                if (Constants.ENABLE_DEBUG &&
                    kdTrees[firstNonFullTree].insert(coordinate) == false) {
                        throw new RuntimeException("\n\nfailure to insert coordinate" +
                                                   " in KdTreeLogarithmic.insert\n");
                }
                if (Constants.ENABLE_DEBUG &&
                        kdTrees[firstNonFullTree].insertedNode == null) {
                        throw new RuntimeException("\n\nnull inserted node" +
                                                   " in KdTreeLogarithmic.insert\n");
                }

                // Create links between the new AVL node and the new k-d node,
                // and return because there is no need to build a new k-d tree.
                avlTree.insertedNode.kdTreeIndex = (short) firstNonFullTree;
                avlTree.insertedNode.kdTreeNode = kdTrees[firstNonFullTree].insertedNode;
                kdTrees[firstNonFullTree].insertedNode.avlTreeNode = avlTree.insertedNode;

                // Verify that the k-d tree contains the correct number of nodes; 
                if ( Constants.ENABLE_DEBUG &&
                     (kdTrees[firstNonFullTree].size() > (1 << firstNonFullTree) ) ) {
                        throw new RuntimeException("\n\ntree " + firstNonFullTree +
                                                   " size = " + kdTrees[firstNonFullTree].size() +
                                                   " in simple insert of KdTreeLogarithmic.insert\n");
                }

                // Walk the list of k-d tree nodes and ensure that the corresponding
                // AVL tree node indexes the newly created k-d tree in the kdTrees array.
                KdNode node = kdTrees[firstNonFullTree].getHead();
                while (node != null) {
                    node.avlTreeNode.kdTreeIndex = (short) firstNonFullTree;
                 }
            }
            return true;
        }

        // Either insertion of the new key-value pair into a non-empty tree
        // was not requested or a non-empty tree was not found.
        //
        // Scan the k-d trees, beginning with B_1, find the smallest tree
        // that is empty. Calculate the sum of tree sizes |B_i|, including
        // |B_0|. Do not allow B_0 to be identified as the smallest empty
        // tree B_i because B_(i-1) might be used instead of B_i below.
        int firstEmptyTree = -1;
        int powerOf2 = 1; // Initialized to 2^(i-1) for i loop index below.
        int treeSizeSum = treeSize(kdTrees[0]);
        for (int i = 1; i < Constants.MAX_POWER_OF_2; ++i) {
            if (treeSize(kdTrees[i]) == 0) {
                firstEmptyTree = i;
                break;
            }
            treeSizeSum += treeSize(kdTrees[i]);
            powerOf2 <<= 1;
        }

        // If no empty tree was found, signal an error.
        if (firstEmptyTree < 0) {
            throw new RuntimeException("\n\nfailed to find empty tree" +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // The empty tree |B_i| was found, so perform the following comparison.
        //
        //     |B_0| + |B_1| + ... + |B_(i-1)| + 1 > 2^(i-1)
        //
        // If the above comparison is true, build a B_i from the key-value
        // pair plus the contents of all trees up to and including B_(i-1).
        // If the above comparison is false, build a B_(i-1) instead.
        // TreeSizeSum is incremented to count the inserted key-value pair.
        int treeIndex = (++treeSizeSum > powerOf2) ? firstEmptyTree : firstEmptyTree - 1;

        // Create a new k-d tree with a single k-d node for the inserted key-value pair.
        KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                               executor,
                                               maxSubmitDepth);

        if (Constants.ENABLE_DEBUG && tree.insert(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to insert coordinate into" +
                                       " k-d tree in KdTreeLogarithmic.insert\n");
        }
        if (Constants.ENABLE_DEBUG && tree.insertedNode == null) {
            throw new RuntimeException("\n\ninserted node is null" +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // Create links between the new AVL node and the single k-d node.
        avlTree.insertedNode.kdTreeNode = tree.insertedNode;
        tree.insertedNode.avlTreeNode = avlTree.insertedNode;

        // Append the k-d nodes from trees B_0 ... B_(i-1) to the node list of
        // the newly created k-d tree. Destroy trees B_0 ... B_(i-1) by setting
        // their kdTrees[] elements to null. B_0 may already be null; however,
        // B_i for i > 0 cannot be null although it may be empty.
        if (kdTrees[0] != null) {
            tree.addList(kdTrees[0]);
            kdTrees[0] = null;
        }
        for (int i = 1; i < firstEmptyTree; ++i) {
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
        // treeSizeSum was incremented above to include the inserted key-value pair.
        if (Constants.ENABLE_DEBUG && tree.size() != treeSizeSum) {
            throw new RuntimeException("\n\ntree size = " + kdTrees[treeIndex].size() +
                                       "  != list size = " + treeSizeSum +
                                       " in KdTreeLogarithmic.insert\n");
        }
        if ( Constants.ENABLE_DEBUG &&
             ( (tree.size() <= (1 << (firstEmptyTree - 2) ) ) || 
               (tree.size() > (1 << firstEmptyTree) ) ) ) {
                   throw new RuntimeException("\n\ntree " + firstEmptyTree +
                                              " size = " + tree.size() + "in elaborate" +
                                              " insert of KdTreeLogarithmic.insert\n");
        }

        // Walk the list of k-d tree nodes and ensure that the corresponding
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
        if (Constants.ENABLE_DEBUG &&
            kdTrees[kdTreeIndex].deletedNode == null) {
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
        // In the parlance of that discussion, the key-value pair "p" has
        // already been "deleted" from B_i (aka B_k), and the "dictionary"
        // (i.e., the AVL tree) has already been updated, above.
        //
        // A departure from the above-reference discussion by Overmars
        // and van Leeuwen is that deletion of the key-value pair has
        // been followed by rebalancing B_i, as suggested in their article,
        // "Dynamic multi-dimensional data structures based on quad- and
        // k-d trees," Acta Informatica 17: 267-285, 1982.
        //
        // According to Overmars and van Leeuwen, no further processing is
        // necessary if, after deletion of the key-value pair, the following
        // post-deletion constraint is satisfied.
        // 
        //     B_i is empty or |B_i| > |2^(i-2)|
        //
        // This post-deletion constraint ignores the constraint |B_1| <= 2^i
        // because the tree-insertion and tree-deletion algorithms proposed
        // by Overmars and van Leeuwen guarantee that this constraint is
        // satisfied. However, some of the examples discussed immediately
        // below use the requirement that |B_i| <= 2^i in their logic.
        //
        // For non-empty B_2, the post-deletion constraint requires that
        // |B_2| > 2^0 and post-deletion, |B_2| can be 3, 2, 1, or 0.
        //
        // For non-empty B_0 before deletion, |B_0| = 2^0 so deletion
        // of the key-value pair produces an empty tree that satisfies
        // the post-deletion constraint.
        //
        // For non-empty B_1, the post-deletion constraint requires that
        // |B_1| > 2^(-1) and it is unclear how to interpret 2^(-1) in
        // the context of Overmars' and van Leeuwen's constraint. However,
        // in general, 2^(i-1) > 2^(i-2), so |B_i| = 2^(i-1) > 2^(i-2)
        // and hence Overmars' and van Leeuwen's constraint is satisfied
        // if |B_i| = 2^(i-1). After deletion, B_1 is either empty or it
        // contains one key-value pair. Hence, |B_1| = 2^0 > 2^(-1) so
        // non-empty, post-deletion B_1 satifies the post-deletion constraint.
        //
        // The above examples demonstrate that it is unnecessary to test
        // whether B_0 and B_1 satisfy the post-deletion constraint
        // |B_i| > 2^(i-2). But it is necessary to test that constraint
        // for B_i where i >= 2.
        if (kdTreeIndex == 0 || kdTreeIndex == 1 ||
            kdTrees[kdTreeIndex].size() > (1 << (kdTreeIndex - 2) ) )
        {
            return true;
        }

        // |B_i| <= 2^(i-2) and in fact, for a correctly built tree,
        // |B_i| = 2^(i-2) after deletion, because before deletion
        // |B_i| > 2^(i-2) so create a new, empty k-d tree for building.
        KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                               executor,
                                               maxSubmitDepth);

        // Overmars' and van Leeuwen's algorithm determines which
        // k-d tree to build and from which trees to cull the nodes
        // used to build the tree. Destroy the trees from which the
        // nodes are culled by setting their elements in the kdTrees
        // array to null.
        int treeIndex;
        // Is B_(i-1) empty?
        if (treeSize(kdTrees[kdTreeIndex-1]) > 0) {
            // B_(i-1) is not empty, so the k-d nodes are
            // culled from B_(i-1) and B_i
            tree.addList(kdTrees[kdTreeIndex-1]);
            tree.addList(kdTrees[kdTreeIndex]);
            kdTrees[kdTreeIndex-1] = kdTrees[kdTreeIndex] = null;
            // Is |B_(i-1) > 2^(i-2)?
            if (kdTrees[kdTreeIndex-1].size() > (1 << (kdTreeIndex-2))) {
                // |B_(i-1)| > 2^(i-2) and, as discussed above, |B_i| = 2^(i-2)
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
        } else if (treeSize(kdTrees[kdTreeIndex-2]) > 0) {
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
            // B_(i-2) using the k-d nodes culled from B_i
            treeIndex = kdTreeIndex - 2;
            tree.addList(kdTrees[kdTreeIndex]);
            kdTrees[kdTreeIndex] = null;
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

        // Verify that the k-d tree contains the correct number of nodes.
        if ( Constants.ENABLE_DEBUG &&
             ( (tree.size() <= (1 << (treeIndex - 2) ) ) || 
               (tree.size() > (1 << treeIndex) ) ) ) {
                   throw new RuntimeException("\n\ntree " + treeIndex +
                                              " size = " + tree.size() +
                                              " in KdTreeLogarithmic.erase\n");
        }

        // Walk the list of k-d tree nodes and ensure that the corresponding
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
            if (treeSize(kdTrees[i]) > 0) {
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
                count += super.verifyKdTree();
            }
        }

        return count;
    }
} // class KdTreeLogarithmic
