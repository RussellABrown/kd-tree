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

        // Create histograms to count rebalancing operations;
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

        // Insertion did not create a new AVL node, so the key was a duplicate;
        // hence, the value was added to a k-d node's values set in a k-d tree,
        // and no further processing is necessary.
        if (avlTree.insertedNode != null) {
            return true;
        }

        // Insertion created a new AVL node, so the key-value pair will be
        // inserted into a k-d tree in the manner discussed by Overmars and
        // van Leeuwen in "Two general methods for dynamizing decomposable
        // searching problems," Computing 26:155-166, 1981. In particular,
        // the size |B_i| of the ith k-d tree B_i is constrained to lie
        // in the range 2^(i-2) < |B_i| <= 2^i, unless the tree is empty
        // in which case the size is zero.
        //
        // In a departure from the above-referenced discussion by Overmars
        // and van Leeuwen, do not perforce build a new k-d tree to insert
        // the key-value pair; instead, insert into a non-empty k-d tree if
        // that tree statisfies the condition |B_i| < 2^i
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
               kdTrees[firstNonFullTree].insert(coordinate);

               // Verify that the key-value pair was inserted.
               if (kdTrees[firstNonFullTree].insertedNode == null) {
                  throw new RuntimeException("\n\nfailure to insert coordinate" +
                                             " in KdTreeLogarithmic.insert\n");
               }

                // Create links between the new AVL node and the new k-d node, and return.
                avlTree.insertedNode.kdTreeIndex = (short) firstNonFullTree;
                avlTree.insertedNode.kdTreeNode = kdTrees[firstNonFullTree].insertedNode;
                kdTrees[firstNonFullTree].insertedNode.avlTreeNode = avlTree.insertedNode;

                return true;
           }
       }

        // Either insertion of the new key-value pair into a non-empty tree
        // was not requested or a non-empty tree was not found.
        //
        // Scan the k-d trees, beginning with the second smallest tree, to
        // find the smallest tree that is empty. Calculate the sum of tree
        // sizes |B_i|, including the size of the smallest tree.
        int firstEmptyTree = -1;
        int powerOf2 = 1; // Initialized to 2^(i-1) for i computed below.
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
        // NOTE: treeSizeSum is incremented to include the key-value pair.
        int treeIndex = (++treeSizeSum > powerOf2) ? firstEmptyTree : firstEmptyTree - 1;

        // Create a k-d tree for the key-value pair, and append the k-d nodes
        // from trees B_0 ... B_(i-1) to the node list of that tree. Empty
        // trees B_0 ... B_(i-1) by setting their kdTrees[] elements to null.
        KdTreeDynamic tree = KdTreeDynamic.createKdTree(coordinate,
                                                        executor,
                                                        maxSubmitDepth,
                                                        insertionHistogramLog);
        for (int i = 0; i < firstEmptyTree; ++i) {
            tree.addList(kdTrees[i]);
            kdTrees[i] = null;
        }

        // Create links between the new AVL node and the k-d node that
        // contains the key-value pair.
        avlTree.insertedNode.kdTreeNode = tree.root;
        tree.root.avlTreeNode = avlTree.insertedNode;

        // Create a k-d tree from the list of k-d nodes.
        kdTrees[treeIndex] = KdTreeDynamic.createKdTree(tree,
                                                        executor,
                                                        maxSubmitDepth,
                                                        insertionHistogramLog);

        // Verify that the k-d tree includes the correct number of nodes.
        // NOTE: treeSizeSum was incremented above to include the key-value pair.
        if (kdTrees[treeIndex].size() != treeSizeSum ) {
            throw new RuntimeException("\n\ntree size = " + kdTrees[treeIndex].size() +
                                       "  != list size = " + treeSizeSum +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // Walk the list of k-d tree nodes and ensure that the corresponding
        // AVL tree node indexes the newly created k-d tree.
        KdNode node = kdTrees[treeIndex].getHead();
        while (node != null) {
            node.avlTreeNode.kdTreeIndex = (short) treeIndex;
        }

        return true;
    }

} // class KdTreeLogarithmic
