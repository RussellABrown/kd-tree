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
        if (coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null in KdTreeLogarithmic.insert\n");
        }

        // Insert the key-value pair into the AVL Tree.
        if (avlTree.insert( coordinate.getKey(), coordinate.getValue() ) == false) {
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
        // Scan the k-d trees, beginning with the smallest tree, to find
        // the smallest tree for which |B_i| < 2^i
        if (Constants.ENABLE_NON_FULL_INSERTION) {
            int firstNonFullTree = -1;
            int powerOf2 = 1; // Initialized to 2^i for i computed below
            for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
                if (kdTrees[i]!= null && kdTrees[i].size() < powerOf2) {
                    firstNonFullTree = i;
                    break;
                }
                powerOf2 <<= 1;
            }
            // Was a non-full tree fount?
            if (firstNonFullTree >= 0) {
                // Yes, so insert the key-value pair into that tree.
                if (kdTrees[firstNonFullTree].insert(coordinate) == false) {
                    throw new RuntimeException("\n\nfailure to insert coordinate" +
                                               " in KdTreeLogarithmic.insert\n");
                }
                if (kdTrees[firstNonFullTree].insertedNode == null) {
                    throw new RuntimeException("\n\nnull inserted node" +
                                                " in KdTreeLogarithmic.insert\n");
                }

                // Create links between the new AVL node and the new k-d node,
                // and return because there is no need to build a new k-d tree.
                avlTree.insertedNode.kdTreeIndex = (short) firstNonFullTree;
                avlTree.insertedNode.kdTreeNode = kdTrees[firstNonFullTree].insertedNode;
                kdTrees[firstNonFullTree].insertedNode.avlTreeNode = avlTree.insertedNode;

                return true;
           }
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

        if (tree.insert(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to insert coordinate into" +
                                       " k-d tree in KdTreeLogarithmic.insert\n");
        }

        // Create links between the new AVL node and the single k-d node.
        avlTree.insertedNode.kdTreeNode = tree.root;
        tree.root.avlTreeNode = avlTree.insertedNode;

        // Append the k-d nodes from trees B_0 ... B_(i-1) to the node list of
        // the newly created k-d tree. Destroy trees B_0 ... B_(i-1) by setting
        // their kdTrees[] elements to null. B_0 may already be null; however,
        // B_i for i > 0 cannot be null.
        if (kdTrees[0] != null) {
            tree.addList(kdTrees[0]);
            kdTrees[0] = null;
        }
        for (int i = 1; i < firstEmptyTree; ++i) {
            tree.addList(kdTrees[i]);
            kdTrees[i] = null;
        }

        // Create a k-d tree from the list of k-d nodes.
        if (Constants.ENABLE_1TO3 && tree.listSize() <= 3) {
            kdTrees[treeIndex] = KdTreeDynamic.createKdTree1to3(tree,
                                                                executor,
                                                                maxSubmitDepth,
                                                                insertionHistogramLog);
        } else {
            kdTrees[treeIndex] = KdTreeDynamic.createKdTree(tree,
                                                            executor,
                                                            maxSubmitDepth,
                                                            insertionHistogramLog);
        }

        // Verify that the k-d tree includes the correct number of nodes;
        // treeSizeSum was incremented above to include the inserted key-value pair.
        if (kdTrees[treeIndex].size() != treeSizeSum ) {
            throw new RuntimeException("\n\ntree size = " + kdTrees[treeIndex].size() +
                                       "  != list size = " + treeSizeSum +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // Walk the list of k-d tree nodes and ensure that the corresponding
        // AVL tree node indexes the newly created k-d tree in the kdTrees array.
        KdNode node = kdTrees[treeIndex].getHead();
        while (node != null) {
            node.avlTreeNode.kdTreeIndex = (short) treeIndex;
        }

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
        if (coordinate == null) {
            throw new RuntimeException("\n\ncoordinate is null in KdTreeLogarithmic.erase\n");
        }

        // Attempt to erase the key-value pair from the AVL tree,
        // and return false if the key was not in the tree.
        if (avlTree.erase( coordinate.getKey(), coordinate.getValue() ) == false) {
            return false;
        }

        // If erasure did not remove an AVL node, the key mapped to more than
        // one value; hence, a value was removed from a k-d node's values set
        // but no k-d node was removed from its k-d tree. No further processing
        // is necessary.
        if (avlTree.removedNode == null) {
            return true;
        }

        // Erasure removed an AVL node, so remove the corresponding k-d node
        // from its k-d tree.
        if (kdTrees[avlTree.removedNode.kdTreeIndex].erase(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to erase coordinate from" +
                                       " k-d tree in KdTreeLogarithmic.erase\n");
        }
        if (kdTrees[avlTree.removedNode.kdTreeIndex].deletedNode == null) {
            throw new RuntimeException("\n\nnull deleted node" +
                                        " in KdTreeLogarithmic.erase\n");
        }

        // Verify that the removed AVL and k-d nodes are correctly linked.
        if ( (avlTree.removedNode.kdTreeNode !=
              kdTrees[avlTree.removedNode.kdTreeIndex].deletedNode) ||
             (avlTree.removedNode !=
              kdTrees[avlTree.removedNode.kdTreeIndex].deletedNode.avlTreeNode) )
        {
            throw new RuntimeException("\n\nAVL and k-d nodes incorrectly" +
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
        // already been "deleted" from B_i (aka B_k) and the "dictionary"
        // (i.e., the AVL tree) has already been updated.
        //
        // A departure from the above-reference discussion by Overmars
        // and van Leeuwen is that "deletion" of "p" was followed by
        // rebalancing of "B_k", which was not described and apparently
        // not contemplated by Overmars and van Leeuwen, although they
        // suggested rebalancing in their subsequent article, "Dynamic
        // multi-dimensional data structures based on quad- and k-d trees,"
        // Acta Informatica 17: 267-285, 1982.
        //
        // According to Overmars and van Leeuwen, no further processing is
        // necessary if, after "deletion" of "p," the following post-deletion
        // condition is satisfied.
        // 
        //     B_i is empty or |2^(i-2)| < |B^i|
        //
        // This post-deletion condition ignores the constraint |B_1| <= 2^i
        // because the tree-insertion and tree-deletion algorithms proposed
        // by Overmars and van Leeuwen guarantee that this constraint is
        // satisfied. However, some of the examples discussed immediately
        // below use the requirement that |B_i| <= 2^i in their logic.
        //
        // For non-empty B_2, the post-deletion condition requires that
        // |B_2| > 2^0.
        //
        // For non-empty B_0, |B_0| <= 2^0 so deletion of the one allowed
        // k-d node produces an empty tree that satisfies the post-deletion
        // condition.
        //
        // For non-empty B_1, the post-deletion condition requires that
        // |B_1| > 2^(-1) and it is unclear how to interpret 2^(-1) in
        // the context of Overmars' and van Leeuwen's proposed constraints.
        // However, because 2^(i-1) > 2^(i-2), |B_i| = 2^(i-1) > 2^(i-2),
        // and further because after deletion, B_1 is either empty or it
        // contains one k-d node. Hence, |B_1| = 2^0 > 2^(-1), and so
        // B_1 satifies the post-deletion condition.
        //
        // The above examples demonstrate that it is unnecessary to test
        // whether B_0 and B_1 satisfy the post-deletion condition. But
        // necessary to test that condition for B_i where i >= 2.
        if (avlTree.removedNode.kdTreeIndex == 0 ||
            avlTree.removedNode.kdTreeIndex == 1 ||
            kdTrees[avlTree.removedNode.kdTreeIndex].size() >
                   (1 << (avlTree.removedNode.kdTreeIndex - 2) ) )
        {
            return true;
        }

        return true;
    }

} // class KdTreeLogarithmic
