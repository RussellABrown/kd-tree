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

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
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
    protected long[] insertHistogramDyn, deleteHistogramDyn, swapHistogramLog;

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
        insertionHistogramLog = new long[Constants.MAX_POWER_OF_2];
        deletionHistogramLog = new long[Constants.MAX_POWER_OF_2];
        swapHistogramLog = new long[Constants.MAX_POWER_OF_2];
        insertHistogramDyn = new long[Constants.MAX_POWER_OF_2];
        deleteHistogramDyn = new long[Constants.MAX_POWER_OF_2];

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
        if (!inserted) {
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
        // unless the tree is empty, where the allowed size |B_i| = 0.
        //
        // The above equation constrains tree sizes as follows, where
        // for B_0 and B_1, 2^(i-2) is undefined so a lower limit of
        // 0 instead of 2^(i-2) is chosen for these trees.
        //
        // |B_3| = 0 or 2 < |B_3| <= 8
        // |B_2| = 0 or 1 < |B_2| <= 4
        // 0 <= |B_1| <= 2
        // 0 <= |B_0| <= 1
        // 
        // In a departure from the above-referenced discussion by Overmars
        // and van Leeuwen, it is not always necessary to build a new tree
        // to insert the key-value pair; instead, that pair may be inserted
        // into a non-full (aka sparse) tree that statisfies the constraints
        // 
        //     2^(i-2) < 1 + |B_i| <= 2^i
        //
        // where 1 is the size of the key-value pair and |B_i| is the size
        // of the sparse tree BEFORE insertion of the key-value pair.
        // Subtracting 1 from all three terms of this equation yields
        //
        //     2^(i-2) <= |B_i| < 2^i
        //
        // However, prior to insertion, the constraint 2^(i-2) < |B_i|
        // is imposed by the insertion and erasure algorithms, and
        // limits the pre-insertion, sparse tree sizes as follows.
        //
        //     2^(i-2) < |B_i| < 2^i
        //
        // The above equation constrains the pre-insertion treee sizes
        // as follows.
        //
        // 2 < |B_3| < 8
        // 1 < |B_2| < 4
        // 0 <= |B_1| < 2
        // 0 <= |B_0| < 1
        //
        // Attempt to insert the key-value pair into B_0 or B_1.
        int insertTree = -1;
        if (getSize(kdTrees[0]) == 0) {
            insertTree = 0;
        } else if (getSize(kdTrees[1]) == 0 || getSize(kdTrees[1]) == 1) {
            insertTree = 1;
        } else if (Constants.ENABLE_DEBUG &&
                   (getSize(kdTrees[0]) < 0 || getSize(kdTrees[0]) > 1))
        {
            throw new RuntimeException("\n\ntree 0 " + " size = " +
                                        getSize(kdTrees[0]) +
                                        " in KdTreeLogarithmic.insert\n");
        } else if (Constants.ENABLE_DEBUG &&
                   (getSize(kdTrees[1]) < 0 || getSize(kdTrees[1]) > 2))
        {
            throw new RuntimeException("\n\ntree 1 " + " size = " +
                                        getSize(kdTrees[1]) +
                                        " in KdTreeLogarithmic.insert\n");
        }

        // If possible, insert the key-value pair into B_0 or B_1?
        if (insertTree >= 0) {
            // Yes, so insert the key-value pair into the tree, which may be null
            // because the getSize method returns 0 for a null or empty tree.
            if (kdTrees[insertTree] == null) {
                kdTrees[insertTree] = new KdTreeDynamic(coordinate.getKey().length,
                                                        executor,
                                                        maxSubmitDepth);
            }
            // Insertion into the tree, and optionally rebalancing that tree, are
            // performed as described by Brown in "A dynamic, self-balancing
            // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
            inserted = kdTrees[insertTree].insert(coordinate);
            if (!inserted) {
                throw new RuntimeException("\n\nfailure to simply insert" +
                                            " coordinate into tree" + insertTree +
                                            " in KdTreeLogarithmic.insert\n");
            }
            if (kdTrees[insertTree].insertedNode == null) {
                throw new RuntimeException("\n\nnull simply inserted k-d node " +
                                            " in tree" + insertTree +
                                            " of KdTreeLogarithmic.insert\n");
            }

            // Create links between the new AVL node and the new k-d node,
            // set the AVL tree's k-d tree index to the k-d tree, and
            // return because there is no need to build a new k-d tree.
            avlTree.insertedNode.kdTreeNode = kdTrees[insertTree].insertedNode;
            kdTrees[insertTree].insertedNode.avlTreeNode = avlTree.insertedNode;
            avlTree.insertedNode.kdTreeIndex = (short) insertTree;

            // Verify that the size of the tree does not exceed its limits.
            if (Constants.ENABLE_DEBUG &&
                ((insertTree == 0 &&
                  (getSize(kdTrees[0]) <= 0 || getSize(kdTrees[0]) > 1)) ||
                 (insertTree == 1 &&
                  (getSize(kdTrees[1]) <= 0 || getSize(kdTrees[1]) > 2))))
            {
                throw new RuntimeException("\n\ntree " + insertTree +" size = " +
                                            getSize(kdTrees[insertTree]) + " after simple" +
                                            " insert of KdTreeLogarithmic.insert\n");
            }
            return true;
        }

        // Insertion failed for B_0 and B_1, so inspect B_i for i >= 2
        // in an attempt to find the smallest sparse tree that satisfies
        // the constraint 2^(i-2) < |B_i| < 2^i
        //
        // NOTE: it appears that the search for a sparse tree ought to
        // preceed the search (below) for an empty tree; otherwise, a
        // sparse tree might not be found.
        if (Constants.ENABLE_SPARSE_INSERTION) {
            int sparseTree = -1;
            // Initialize power-of-2 limits for i = 2 and exclude B_0 and B_1
            // from the following 'for' loop.
            int powerOf2Lower = 1; // Initialized to 2^(i-2) for i loop index below.
            int powerOf2Upper = 4; // Initialized to 2^i for i loop index below.
            for (int i = 2; i < Constants.MAX_POWER_OF_2; ++i) {
                if (getSize(kdTrees[i]) > powerOf2Lower &&
                    getSize(kdTrees[i]) < powerOf2Upper) // No need to check powerOf2Upper
                {
                    sparseTree = i;
                    break;
                } else if (Constants.ENABLE_DEBUG && sparseTree >= 0 &&
                           (getSize(kdTrees[sparseTree]) <= powerOf2Lower ||
                            getSize(kdTrees[sparseTree]) > powerOf2Upper))
                { 
                    throw new RuntimeException("\n\ntree " + sparseTree + " size = " +
                                                getSize(kdTrees[sparseTree]) +
                                                " before sparse" +
                                                " insert of KdTreeLogarithmic.insert\n");
                }
                powerOf2Lower <<= 1;
                powerOf2Upper <<= 1;
            }
 
            // Was a sparse tree found wherein 2^(i-2) < |B_i| < 2^i ? 
            if (sparseTree >= 0) {
                // Yes, so insert the key-value pair into that tree, which may be null.
                if (kdTrees[sparseTree] == null) {
                    kdTrees[sparseTree] = new KdTreeDynamic(coordinate.getKey().length,
                                                            executor,
                                                            maxSubmitDepth);
                }
                // Insertion into the sparse tree, and optionally rebalancing that tree,
                // are performed as described by Brown in "A dynamic, self-balancing
                // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
                inserted = kdTrees[sparseTree].insert(coordinate);
                if (!inserted) {
                    throw new RuntimeException("\n\nfailure to insert coordinate" +
                                                " into sparse tree" + sparseTree +
                                                " in KdTreeLogarithmic.insert\n");
                }
                if (kdTrees[sparseTree].insertedNode == null) {
                    throw new RuntimeException("\n\nnull inserted k-d node " +
                                               " in sparse tree" + sparseTree +
                                               " of KdTreeLogarithmic.insert\n");
                }

                // Create links between the new AVL node and the new k-d node,
                // set the AVL tree's k-d tree index to the sparse tree, and
                // return because there is no need to build a new k-d tree.
                avlTree.insertedNode.kdTreeNode = kdTrees[sparseTree].insertedNode;
                kdTrees[sparseTree].insertedNode.avlTreeNode = avlTree.insertedNode;
                avlTree.insertedNode.kdTreeIndex = (short) sparseTree;

                // Verify that the size of the tree does not exceed its limits.
                if (Constants.ENABLE_DEBUG &&
                    (getSize(kdTrees[sparseTree]) <= (1 << (sparseTree-2)) ||
                     getSize(kdTrees[sparseTree]) > (1 << sparseTree)))
                {
                    throw new RuntimeException("\n\ntree " + sparseTree +" size = " +
                                               getSize(kdTrees[sparseTree]) + " after simple" +
                                               " insert of KdTreeLogarithmic.insert\n");
                }
                return true;
            }
        }

        // Insertion of the new key-value pair into a sparse tree wherein 2^(i-2) < |B_i| < 2^i
        // either was not requested or a sparse tree was not found, so find the smallest
        // empty tree, and replace either that tree or the next smaller tree with a new tree.
        //
        // Scan B_i for i >= 2 to find the smallest empty tree. Calculate
        // the sum of tree sizes |B_(i-1)| including |B_0| and |B_1|.
        int emptyTree = -1;
        // treeSizeSum includes the size of the to-be-inserted key-value
        // pair (1) as well as |B_0| and |B_1|.
        int treeSizeSum = 1 + getSize(kdTrees[0]) + getSize(kdTrees[1]);
        int powerOf2 = 2; // Initialized to 2^(i-1) for i loop index below.
        for (int i = 2; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) <= 0) {
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
        int treeIndex = (treeSizeSum > powerOf2) ? emptyTree : emptyTree - 1;

        // Create a new k-d tree and insert the key-value pair.
        final KdTreeDynamic tree = new KdTreeDynamic(coordinate.getKey().length,
                                                     executor,
                                                     maxSubmitDepth);

        inserted = tree.insert(coordinate);
        if (Constants.ENABLE_DEBUG && !inserted) {
            throw new RuntimeException("\n\nfailed to insert coordinate into" +
                                        " nascent tree " + treeIndex +
                                        " in KdTreeLogarithmic.insert\n");
        }
        if (tree.insertedNode == null) {
            throw new RuntimeException("\n\ninserted node is null" +
                                        " in nascent tree " + treeIndex +
                                        " of KdTreeLogarithmic.insert\n");
        }
        final KdNode saveInsertedNode = tree.insertedNode;

        // Append the k-d nodes from trees B_0 ... B_(i-1) to the node list of
        // the newly created k-d tree, and then clear trees B_0 ... B_(i-1)
        // so that they appear empty.
        for (int i = 0; i < emptyTree; ++i) {
            tree.addList(kdTrees[i]);
            if (kdTrees[i] != null) {
                kdTrees[i].clear();
            }
        }

        // Replace the new k-d tree with a k-d tree created from the list of k-d nodes,
        // and set its reference to the inserted k-d node from previously created tree
        // into which the key-value pair was inserted, because the createKdTree method
        // does not set the reference to this inserted node.
        kdTrees[treeIndex] = createKdTree(tree,
                                          numDimensions,
                                          executor,
                                          maxSubmitDepth,
                                          insertionHistogramLog);

        kdTrees[treeIndex].insertedNode = saveInsertedNode;

        // Verify that the k-d tree contains the correct number of nodes; 
        if (Constants.ENABLE_DEBUG &&
            getSize(kdTrees[treeIndex]) != treeSizeSum)
        {
            throw new RuntimeException("\n\ntree size = " + getSize(kdTrees[treeIndex]) +
                                       "  !=  list size = " + treeSizeSum +
                                       " in KdTreeLogarithmic.insert\n");
        }

        // Verify the correct size of the k-d tree.
        if (Constants.ENABLE_DEBUG &&
            (getSize(kdTrees[treeIndex]) <= (1 << (treeIndex-2))  ||
             getSize(kdTrees[treeIndex]) > (1 << treeIndex)))
        {
            throw new RuntimeException("\n\ntree " + treeIndex + " size = " +
                                       getSize(kdTrees[treeIndex]) + "in elaborate insert" +
                                       " of KdTreeLogarithmic.insert\n");
        }

        // Create links between the new AVL node and the inserted k-d node.
        avlTree.insertedNode.kdTreeNode = kdTrees[treeIndex].insertedNode;
        kdTrees[treeIndex].insertedNode.avlTreeNode = avlTree.insertedNode;

        // Walk the list of k-d tree nodes and ensure that each corresponding
        // AVL tree node indexes the correct k-d tree in the kdTrees array.
        KdNode node = kdTrees[treeIndex].head;
        while (node != null) {
            node.avlTreeNode.kdTreeIndex = (short) treeIndex;
            node = node.next;
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
             (kdTreeIndex < 0) || (kdTreeIndex >= Constants.MAX_POWER_OF_2) )
        {
            throw new RuntimeException("\n\nk-d tree index = " + kdTreeIndex +
                                       " is out of range in KdTreeLogarithmic.erase\n");
        }
        if (kdTrees[kdTreeIndex].erase(coordinate) == false) {
            throw new RuntimeException("\n\nfailed to erase coordinate from" +
                                       " k-d tree in KdTreeLogarithmic.erase\n");
        }
        if (kdTrees[kdTreeIndex].deletedNode == null) {
            throw new RuntimeException("\n\nnull deleted node" +
                                       " in KdTreeLogarithmic.erase\n");
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
        // unless the tree is empty, where the allowed size |B_i| = 0.
        //
        // The above equation constrains tree sizes as follows, where
        // for B_0 and B_1, 2^(i-2) is undefined so a lower limit of
        // 0 instead of 2^(i-2) is chosen for these trees.
        // 
        // |B_3| = 0 or 2 < |B_3| <= 8
        // |B_2| = 0 or 1 < |B_2| <= 4
        // 0 <= |B_1| <= 2
        // 0 <= |B_0| <= 1
        // 
        // After removal of a node, the tree size is constrained as follows.
        // 
        //     2^(i-2) < 1 + |B_i| <= 2^i
        //
        // where 1 is the size of the key-value pair and |B_i| is the size
        // of the tree AFTER removal of the key-value pair. Subtracting 1
        // from all three terms of this equation yields the following
        // constraints.
        //
        //     2^(i-2) <= |B_i| < 2^i
        //
        // However, for B_0 and B_1, 2^(i-2) is undefined so a lower
        // limit of 0 instead of 2^(i-2) is chosen for these trees 
        // that limits post-removal tree sizes as follows.
        //
        // 2 <= |B_3| < 8
        // 1 <= |B_2| < 4
        // 0 <= |B_1| < 2
        // 0 <= |B_0| < 1
        //
        // A departure from the above-reference discussion by Overmars
        // and van Leeuwen is that deletion of the key-value pair may
        // optionally be followed by rebalancing B_i, as suggested in
        // their article, "Dynamic multi-dimensional data structures based
        // on quad- and k-d trees," Acta Informatica 17: 267-285, 1982.
        //
        // Deletion from the tree, and rebalancing that tree, are
        // performed as described by Brown in "A dynamic, self-balancing
        // k-d tree," www.arxiv.org/abs/2509.08148, 1-19, 2026.
        //
        // According to Overmars and van Leeuwen, no further processing is
        // necessary if, after deletion of the key-value pair, the following
        // post-deletion constraint is satisfied.
        // 
        //     B_i is empty or 2^(i-2) < |B_i|
        //
        // This constraint omits the constraint |B_i| <= 2^i because the
        // tree-insertion algorithms of the insert method guarantee that
        // this contraint is satisfied. However, including this constraint
        // in the algorithm below for B_i where i >= 2, and also consistent
        // with the algorithms of the insert method, the following code
        // handles the cases of B_0 and B_1.
        if (kdTreeIndex == 0) {
            // B_0 should be empty due to erasure of its one allowed node,
            // so no further processing is necessary because |B_0| = 0
            // is allowed by its constraints. Check that B_0 is indeed empty.
            if (Constants.ENABLE_DEBUG && getSize(kdTrees[0]) != 0)
            {
                throw new RuntimeException("\n\ntree 0 size = " +
                                            getSize(kdTrees[0]) +
                                            " after erasure in" +
                                            " KdTreeLogarithmic.erase\n");
            }
            return true;

        } else if (kdTreeIndex == 1) {
            // For a properly built B_1, erasure of one node should not
            // invalidate B_1, so no further processing is necessary
            // because after erasure from a properly built tree,
            // its constraints allow 0 <= |B_1| < 2.
            // Check that B_1 contains the correct number of nodes.
            if (Constants.ENABLE_DEBUG &&
                (getSize(kdTrees[1]) < 0 || getSize(kdTrees[1]) > 1))
            {
                throw new RuntimeException("\n\ntree 1 size = " +
                                            getSize(kdTrees[1]) +
                                            " after erasure in" +
                                            " KdTreeLogarithmic.erase\n");
            }
            return true;

        } else if (getSize(kdTrees[kdTreeIndex]) > (1 << (kdTreeIndex - 2))) {
            // B_i for i >= 2 satisfies the constraint |B_i| > 2^(i-2)
            // so no further processing is necessary.
            return true;

        } else {

            // Given the above processing, |B_i| = 2^(i-2) for i >= 2
            // after erasure of one node because, for a properly built
            // tree, |B_i| > 2^(i-2) before erasure of one node. Hence,
            // failure of that test after erasure of one node (performed
            // immediately above) implies that |B_1| = 2^(i-2).
            // 
            // Declare a tree reference for the tree that will be rebuilt,
            // as well as an index into the kdTrees array for that tree.
            KdTreeDynamic tree;
            int treeIndex;

            // Overmars' and van Leeuwen's algorithm determines which
            // k-d tree to build and from which tree to cull the nodes
            // used to build the tree. After culling, clear the tree
            // from which the nodes are culled so that it is empty.
 
            // Is B_(i-1) empty?
            if (getSize(kdTrees[kdTreeIndex-1]) > 0) {
                // Neither B_(i-1) nor B_i is empty; hence, the nodes
                // are culled from either B_(i-1) or B_i into the new tree.
                //
                // Is |B_(i-1)| > 2^(i-2)?
                if (getSize(kdTrees[kdTreeIndex-1]) > (1 << (kdTreeIndex-2))) {
                    // |B_(i-1)| > 2^(i-2) and |B_i| = 2^(i-2) as discussed
                    // above, so |B_(i-1)| + |B_i| > 2*2^(i-2) = 2^(i-1) and
                    // hence the combined nodes of B_(i-1) and B_i exceed the
                    // the constraint that |B_(i-1)| <= 2^(i-1) so cull the
                    // nodes from B_(i-1) into B_i, and then clear B_(i-1),
                    // because B_i will be rebuilt from the combined nodes
                    // of B_(i-1) and B_i.
                    treeIndex = kdTreeIndex;
                    tree = kdTrees[treeIndex];
                    tree.addList(kdTrees[kdTreeIndex-1]);
                    kdTrees[kdTreeIndex-1].clear();
                } else {
                    // |B_(i-1)| <= 2^(i-2) and |B_i| = 2^(i-2) as discussed
                    // above, so |B_i-1)| + |B_i| <= 2*2^(i-2) = 2^(i-1) and
                    // hence the combined nodes of B_(i-1) and B_i satisfy the
                    // constraint that |B_(i-1)| <= 2^(i-1) so cull the
                    // nodes from B_i into B_(i-1), and then clear B_i,
                    // because B_(i-1) will be rebuilt from the combined
                    // nodes of B_(i-1) and B_i.
                    treeIndex = kdTreeIndex - 1;
                    tree = kdTrees[treeIndex];
                    tree.addList(kdTrees[kdTreeIndex]);
                    kdTrees[kdTreeIndex].clear();
                }
                // Replace the selected tree between B_i and B_(i-1) with the
                // tree created from the combined nodes from B_i and B_(i-1).
                kdTrees[treeIndex] = createKdTree(tree,
                                                  numDimensions,
                                                  executor,
                                                  maxSubmitDepth,
                                                  deletionHistogramLog);

            // Is B_(i-2) empty?
            } else if (getSize(kdTrees[kdTreeIndex-2]) > 0) {
                // Neither B_(i-2) nor B_i is empty but B_(i-1)
                // is empty, so because |B_(i-2)| <= 2^(i-2) and
                // also because |B_i| = 2^(i-2) as discussed above,
                // |B_(i-2)| + |B_i| <= 2*2^(i-2) = 2^(i-1) so
                // cull the nodes from B_i and B_(i-2), clear
                // these trees, and then build B_(i-1) from
                // the combined nodes from B_i and B_(i-2).
                treeIndex = kdTreeIndex - 1;
                tree = kdTrees[treeIndex];
                tree.addList(kdTrees[kdTreeIndex]);
                tree.addList(kdTrees[kdTreeIndex-2]);
                kdTrees[kdTreeIndex].clear();
                kdTrees[kdTreeIndex-2].clear();
                // Replace B_(i-1) with a tree created from
                // the combined nodes from B_i and B_(i-2).
                kdTrees[treeIndex] = createKdTree(tree,
                                                  numDimensions,
                                                  executor,
                                                  maxSubmitDepth,
                                                  deletionHistogramLog);

            } else {
                // B_(i-1) and B_(i-2) are both empty but
                // B_i is not empty and hence, as discussed
                // above, |B_i| = 2^(i-2) so swap B_i and
                // B_(i-2) thus moving the nodes from B_i
                // to B_(i-2) and emptying B_i.
                incrementHistogram(swapHistogramLog, getSize(kdTrees[kdTreeIndex]));
                treeIndex = kdTreeIndex - 2;
                final KdTreeDynamic tmpTree = kdTrees[kdTreeIndex];
                kdTrees[kdTreeIndex] = kdTrees[kdTreeIndex-2];
                kdTrees[kdTreeIndex-2] = tmpTree;
                tree = kdTrees[treeIndex];
            }

            // Verify the correct size of the k-d tree B_i for i >= 2.
            if (Constants.ENABLE_DEBUG && getSize(kdTrees[treeIndex]) >= 2 &&
                (getSize(kdTrees[treeIndex]) <= (1 << (treeIndex-2)) ||
                 getSize(kdTrees[treeIndex]) > (1 << treeIndex)))
            {
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
                node = node.next;
            }

            return true;
        }
    }

    /**
     * <p>
     * The {@code verifyTree} method checks that the children of each node
     * of the k-d tree is correctly sorted relative to that node.
     * </p>
     * 
     * @return the number of nodes in the k-d tree
     */
    protected long verifyKdTree()
    {
        // Verify all k-d trees.
        long count = 0L;
        // Inspect all non-empty k-d trees.
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                // Verify that each non-empty k-d tree B_i conforms
                // to its size constaints 2^(i-2) < |B_i| <= 2^i that
                // require, for example, the following tree sizes.
                //
                // 2 < |B_3| <= 8
                // 1 < |B_2| <= 4
                // 1 <= |B_1| <= 2
                // 0 < |B_0| <= 1
                //
                if ((i == 0 && getSize(kdTrees[i]) > 1) ||
                    (i == 1 && (getSize(kdTrees[i]) < 1 ||
                                getSize(kdTrees[i]) > 2)) ||
                    (i >= 2 &&
                     ((1 << (i-2)) > getSize(kdTrees[i]) ||
                      getSize(kdTrees[i]) > (1 << i))))
                {
                    throw new RuntimeException("\n\ntree " + i + " size = " + 
                                               getSize(kdTrees[i]) +
                                               " in verifyKdTree\n");
                }
                // Verify that each AVL node contains the correct k-d tree index.
                KdNode node = kdTrees[i].head;
                while (node != null) {
                    if (node.avlTreeNode.kdTreeIndex != (short) i) {
                        throw new RuntimeException("\n\nAVL node index = " +
                                                   node.avlTreeNode.kdTreeIndex +
                                                   "  !=  tree index = " + i +
                                                   " in verifyKdTree\n");
                    }
                    node = node.next;
                }
                count += kdTrees[i].verifyKdTree(true);
            }
        }

        return count;
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
        KdNode avlKdTreeNode = null, kdTreeNode = null, kdKdTreeNode = null;
        boolean avlResult = false, kdResult = false;

        // Search the AVL-tree dictionary for the key.
        final AvlNode avlTreeNode = avlTree.find(coordinate);
        if (avlTreeNode != null) {
            // The key was found in an AVL node, so search for the value
            // in the associated k-d node.
            avlKdTreeNode = avlTreeNode.kdTreeNode;
            if (avlKdTreeNode.values.contains(coordinate.getValue())) {
                avlResult = true;
            }
        }

        // Compare the search of the AVL tree to a search of the k-d trees.
        if (Constants.ENABLE_DEBUG) {
            for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
                if (kdTrees[i] != null && kdTrees[i].nodeCount > 0) {
                    kdTreeNode = kdTrees[i].find(coordinate);
                    if (kdTreeNode != null) {
                        kdKdTreeNode = kdTreeNode;
                        // Check whether a prior k-d tree contains the key.
                        if (!kdResult) {
                            kdResult = true;
                        } else {
                            throw new RuntimeException("\n\nkey found twice" +
                                                       " in contains\n");
                        }
                    }
                }
            }
            if (avlResult != kdResult) {
                throw new RuntimeException("\n\nboolean mismatch in contains:" +
                                           " AVL result = " + avlResult +
                                           "  k-d result = " + kdResult + "\n");
            }
            if (avlKdTreeNode != kdKdTreeNode) {
                throw new RuntimeException("\n\nsearch mismatch in contains\n");
            }
        }

        return avlResult;
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
     */
    protected List<KdNode> searchKdTree(final long[] queryLower,
                                        final long[] queryUpper,
                                        final ExecutorService executor,
                                        final int maximumSubmitDepth,
                                        final int p,
                                        final int depth,
                                        final boolean enableAll)
    {
        List<KdNode> nodeList;
        if (Constants.ENABLE_LINKED_LIST) {
            nodeList = new LinkedList<KdNode>();
        } else {
            nodeList = new ArrayList<KdNode>();
        }

        // Search all non-empty k-d trees and append the results.
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                nodeList.addAll(kdTrees[i].searchKdTree(queryLower, queryUpper,
                                                        executor, maximumSubmitDepth,
                                                        p, depth, enableAll));
            }
        }
      return nodeList;
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
     */
    protected List<KdNode> searchKdTree(final long[] queryLower,
                                        final long[] queryUpper,
                                        final ExecutorService executor,
                                        final int maximumSubmitDepth,
                                        final int p,
                                        final int depth,
                                        final boolean[] enable)
    {
        List<KdNode> nodeList;
        if (Constants.ENABLE_LINKED_LIST) {
            nodeList = new LinkedList<KdNode>();
        } else {
            nodeList = new ArrayList<KdNode>();
        }
        
        // Search all non-empty k-d trees and append the results.
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                nodeList.addAll(kdTrees[i].searchKdTree(queryLower, queryUpper,
                                                        executor, maximumSubmitDepth,
                                                        p, depth, enable));
            }
        }
      return nodeList;
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
     * @return a {@link java.util.List List}{@code <}{@link Paire}{@code double, KdNode>>}
     */
    protected List<Paire> findNearestNeighbors(final long[] query,
                                               final int numNeighbors)
    {
        // Search all non-empty k-d trees and add the resulting Paires to the initial list.
        final LinkedList<Paire> initialNeighbors = new LinkedList<Paire>();
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                initialNeighbors.addAll(kdTrees[i].root.findNearestNeighbors(query,
                                                                             numNeighbors));
            }
        }

        // Create a heap and add the Paires from the initial list.
        final NearestNeighborHeap heap = new NearestNeighborHeap(query,
                                                                 initialNeighbors.size());
        for (Paire item : initialNeighbors) {
            heap.add(item.getValue());
        }

        // Append only up to M nearest-neighbor Paires to the final list.
        final LinkedList<Paire> finalNeighbors = new LinkedList<Paire>();
        final int numFewNeighbors =
            (heap.heapDepth() < numNeighbors) ? heap.heapDepth() : numNeighbors;
        for (int i = 0; i < numFewNeighbors; ++i) {
            finalNeighbors.addFirst(heap.removeTop());
        }

        return finalNeighbors;
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
     * @return a {@link java.util.List List}{@code <}{@link Paire}{@code double, KdNode>>}
     */
    protected List<Paire> findNearestNeighbors(final long[] query,
                                               final int numNeighbors,
                                               final boolean[] enable)
    {
        // Search all non-empty k-d trees and add the resulting Paires to the initial list.
        final LinkedList<Paire> initialNeighbors = new LinkedList<Paire>();
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                initialNeighbors.addAll(kdTrees[i].root.findNearestNeighbors(query,
                                                                             numNeighbors));
            }
        }

        // Create a heap and add the Paires from the initial list.
        final NearestNeighborHeap heap = new NearestNeighborHeap(query,
                                                                 initialNeighbors.size(),
                                                                 enable);
        for (Paire item : initialNeighbors) {
            heap.add(item.getValue());
        }

        // Append only up to M nearest-neighbor Paires to the final list.
        final LinkedList<Paire> finalNeighbors = new LinkedList<Paire>();
        final int numFewNeighbors =
            (heap.heapDepth() < numNeighbors) ? heap.heapDepth() : numNeighbors;
        for (int i = 0; i < numFewNeighbors; ++i) {
            finalNeighbors.addFirst(heap.removeTop());
        }

        return finalNeighbors;
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
     * @return a {@link java.util.List List}{@code <}{@link Paire}{@code double, KdNode>>}
     */
    protected  List<Paire> findBruteNeighbors(final long[] query,
                                              final int numNeighbors)
    {
        // Search all non-empty k-d trees and add the resulting Paires to the initial list.
        final LinkedList<Paire> initialNeighbors = new LinkedList<Paire>();
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                initialNeighbors.addAll(kdTrees[i].root.findBruteNeighbors(query,
                                                                           numNeighbors));
            }
        }

        // Create a heap and add the Paires from the initial list.
        final NearestNeighborHeap heap = new NearestNeighborHeap(query,
                                                                 initialNeighbors.size());
        for (Paire item : initialNeighbors) {
            heap.add(item.getValue());
        }

        // Append only up to M nearest-neighbor Paires to the final list.
        final LinkedList<Paire> finalNeighbors = new LinkedList<Paire>();
        final int numFewNeighbors =
            (heap.heapDepth() < numNeighbors) ? heap.heapDepth() : numNeighbors;
        for (int i = 0; i < numFewNeighbors; ++i) {
            finalNeighbors.addFirst(heap.removeTop());
        }

        return finalNeighbors;
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
     * @return a {@link java.util.List List}{@code <}{@link Paire}{@code double, KdNode>>}
     */
    protected List<Paire> findBruteNeighbors(final long[] query,
                                             final int numNeighbors,
                                             final boolean[] enable)
    {
        // Search all non-empty k-d trees and add the resulting Paires to the initial list.
        final LinkedList<Paire> initialNeighbors = new LinkedList<Paire>();
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (getSize(kdTrees[i]) > 0) {
                initialNeighbors.addAll(kdTrees[i].root.findBruteNeighbors(query,
                                                                           numNeighbors,
                                                                           enable));
            }
        }

        // Create a heap and add the Paires from the initial list.
        final NearestNeighborHeap heap = new NearestNeighborHeap(query,
                                                                 initialNeighbors.size());
        for (Paire item : initialNeighbors) {
            heap.add(item.getValue());
        }

        // Append only up to M nearest-neighbor Paires to the final list.
        final LinkedList<Paire> finalNeighbors = new LinkedList<Paire>();
        final int numFewNeighbors =
            (heap.heapDepth() < numNeighbors) ? heap.heapDepth() : numNeighbors;
        for (int i = 0; i < numFewNeighbors; ++i) {
            finalNeighbors.addFirst(heap.removeTop());
        }

        return finalNeighbors;
    }
  
    /**
    * <p>
    * The {@code sumHistograms} method sums the dynamic histograms
    * <p>
    */
    protected void sumHistograms()
    {
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            for (int j = 0; j < Constants.MAX_POWER_OF_2; ++j) {
                if (kdTrees[j] != null) {
                    insertHistogramDyn[i] += kdTrees[j].insertionHistogramDyn[i];
                    deleteHistogramDyn[i] += kdTrees[j].deletionHistogramDyn[i];
                }
            }
        }
    }

    /**
     * <p>
     * The {@code getTreeHeight} method returns the height at the root
     * of the highest {@code KdTreeDynamic} of the {@code KdTreeLogarithmic}.
     * 
     * @return the height of the tallest KdTreeDynamic
     * </p>
     */
    protected int getTreeHeight()
    {
        int height = 0;
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            if (KdTreeDynamic.getSize(kdTrees[i])> 0 && kdTrees[i].root.height > height) {
                height = kdTrees[i].root.height;
            }
        }
        return height;
    }

    /**
     * <p>
     * The {@code treeSize} method returns the sum of sizes of all of the
     * {@code KdTreeDynamic} of the {@code KdTreeLogarithmic}.
     * 
     * @return the tree size if the tree is not null; otherwise, 0
     * </p>
     */
    protected static long getSize(final KdTreeLogarithmic tree)
    {
        long nodeCount = 0L;
        for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
            nodeCount += KdTreeDynamic.getSize(tree.kdTrees[i]);
        }
        return nodeCount;
    }

}
 // class KdTreeLogarithmic
