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
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The {@code AvlNode} class stores a k-dimensional key, as well as
 * references to the "less than" and "greater than" sub-trees, to a
 * KdNode instance, a duplicates count, a balance, and an index
 * into the KdTreeLogarithmic.kdTrees array.
 * </p>
 */
public class AvlNode {
    
    protected long[] tuple;
    private short bal;
    protected short kdTreeIndex;
    private int height; // a pathological tree may be 2^31-1 deep
    protected KdNode kdTreeNode;
    protected AvlNode left, right;

    /**
     * <p>
     * AvlNode constructor
     * 
     * @param key - a {@code long}[]
     * <p>
     */
    public AvlNode(final long[] key)
    {
        bal = 0;
        kdTreeIndex = 0;    // Assign after selection of the k-d tree
        kdTreeNode = null;  // Assign after insertion into the k-d tree.
        left = right = null;

        // Make a copy of the key.
        tuple = new long[key.length];
        System.arraycopy(key, 0, tuple, 0, key.length);
    }

    /**
     * <p>
     * The {@code contains} method searches the AVL tree for the existence
     * of a super-key and returns true if the AVL tree contains it.
     * 
     * @param key - a {@code long}[]
     * @return {@code true} if the tree contains the key; otherwise, {@code false}
     * </p>
     */
    protected boolean contains(final long[] key)
    {
        return ( find(key) != null );
    }

    /**
     * <p>
     * The {@code find} method searches the AVL tree iteratively for the existence of
     * a super-key, and returns a reference to the AvlNode instance that contains it.
     * 
     * NOTE: the super-key does not permute cyclically because an AVL tree is unidimensional.
     * 
     * @param key - a {@code long}[]
     * @return a {@code AvlNode} reference if the tree contains the key; otherwise, {@code null}
     * </p>
     */
    protected AvlNode find(final long[] key)
    {
        AvlNode p = this;
        
        while ( p != null ) { // Iterate; don't use recursion
            if ( MergeSort.superKeyCompare(key, p.tuple, 0) < 0 ) {
                p = p.left; // Follow the left branch
            } else if ( MergeSort.superKeyCompare(key, p.tuple, 0) > 0 ) {
                p = p.right; // Follow the right branch
            } else {
                return p; // Found the key, so return the AvlNode reference
            }
        }
        return null; // Didn't find the key, so return null
    }

    /**
     * <p>
     * The {@code insert} method searches the AVL tree for the existence of
     * a super-key. If the super-key is not found, the super-key is inserted
     * as a new AvlNode and the value is added to the AvlNode's {@code TreeSet),
     * AvlTree.insertedNode is assigned a reference to the new AvlNode, and
     * AvlTree.nodeCount is incremented. But if the super-key is found, the
     * value is added to the {@code TreeSet} of the AvlNode that contains
     * the super-key.
     * 
     * This method is adapted from the {@code search} procedure on pp. 220-221
     * of Nicklaus Wirth's textbook, "Algorithms + Data Structures = Programs"
     * (Prentice-Hall, 1976).
     * 
     * NOTE: the super-key does not permute cyclically because an AVL tree is unidimensional.
     * 
     * @param key - a {@code long}[]
     * @param value - a {@code String}
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    protected AvlNode insert(final long[] key,
                             final String value,
                             final AvlTree tree)
    {
        AvlNode p = this;
        
        if ( MergeSort.superKeyCompare(key, p.tuple, 0) < 0 ) { // Search the left branch?
            if ( p.left != null ) {
                p.left = p.left.insert( key, value, tree );
            } else {
                p.left = tree.insertedNode = new AvlNode( key );
                tree.changed = true;
            }
            if ( tree.changed == true ) {
                switch ( p.bal ) {                  // Left branch has grown higher
                    case 1:                         // Balance restored...
                        p.bal = 0;
                        tree.changed = false;       // ...so height change has been eliminated
                        break;
                    case 0:                         // Tree has become more unbalanced
                        p.bal = -1;
                        break;
                    case -1:		                // Tree must be rebalanced
                        AvlNode p1 = p.left;
                        if ( p1.bal == -1 ) {		// Single LL rotation
                            tree.lli++;
                            p.left = p1.right;
                            p1.right = p;
                            p.bal = 0;
                            p = p1;
                        } else {			        // Double LR rotation
                            tree.lri++;
                            AvlNode p2 = p1.right;
                            p1.right = p2.left;
                            p2.left = p1;
                            p.left = p2.right;
                            p2.right = p;
                            if ( p2.bal == -1 ) {
                                p.bal = 1;
                            } else {
                                p.bal = 0;
                            }
                            if ( p2.bal == 1 ) {
                                p1.bal = -1;
                            } else {
                                p1.bal = 0;
                            }
                            p = p2;
                        }
                        p.bal = 0;
                        tree.changed = false;       // Rebalancing has eliminated height change
                        break;
                    default:
                        {
                            throw new RuntimeException("\n" + "bal = " + p.bal + " is out of range in AvlNode.insert right\n");
                        }
                }
            }
        } else if ( MergeSort.superKeyCompare(key, p.tuple, 0) > 0 ) { // Search the right branch?
            if ( p.right != null ) {
                p.right = p.right.insert( key, value, tree );
            } else {
                p.right = tree.insertedNode = new AvlNode( key );
                tree.changed = true;
            }
            if ( tree.changed == true ) {           // Right branch has grown higher
                switch ( p.bal ) {
                    case -1:                        // Balance restored...
                        p.bal = 0;
                        tree.changed = false;       // ...so height change has been eliminated
                        break;
                    case 0:                         // Tree has become more unbalanced
                        p.bal = 1;
                        break;
                    case 1:                         // Tree must be rebalanced
                        AvlNode p1 = p.right;
                        if ( p1.bal == 1 ) {        // Single RR rotation
                            tree.rri++;
                            p.right = p1.left;
                            p1.left = p;
                            p.bal = 0;
                            p = p1;
                        } else {                    // Double RL rotation
                            tree.rli++;
                            AvlNode p2 = p1.left;
                            p1.left = p2.right;
                            p2.right = p1;
                            p.right = p2.left;
                            p2.left = p;
                            if ( p2.bal == 1 ) {
                                p.bal = -1;
                            } else {
                                p.bal = 0;
                            }
                            if ( p2.bal == -1 ) {
                                p1.bal = 1;
                            } else {
                                p1.bal = 0;
                            }
                            p = p2;
                        }
                        p.bal = 0;
                        tree.changed = false;       // Rebalancing has eliminated height change
                        break;
                    default:
                        {
                            throw new RuntimeException("\n" + "bal = " + p.bal + " is out of range in AvlNode.insert left\n");
                        }
                }
            }
        } else { // The tree already contains the key, so add the value to KdTreeNode.values
            if (p.kdTreeNode != null) {
                p.kdTreeNode.values.add(value);
                tree.insertedNode = null;
            } else {
                throw new RuntimeException("\nkdTreeNode reference is null in AvlNode.insert\n");
            }
            tree.changed = false; // A new AvlNode was not inserted, so the height hasn't changed
        }
        return p; // The root of the rebalanced sub-tree
    }

    /**
     * <p>
     * The {@code balanceLeft} method rebalances following deletion of a left AvlNode.
     * 
     * This method is adapted from the {@code balance1} procedure on pp. 223-224 of Nickaus
     * Nicklaus textbook, "Algorithms + Data Structures = Programs" (Prentice-Hall, 1976).
     * 
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    private AvlNode balanceLeft(final AvlTree tree)
    {
        AvlNode p = this;
        
        switch ( p.bal ) {
            case -1:                    // Balance restored
                p.bal = 0;
                break;
            case 0:                     // Tree has become more unbalanced
                p.bal = 1;
                tree.changed = false;
                break;
            case 1:                     // Tree must be rebalanced
                AvlNode p1 = p.right;
                if ( p1.bal >= 0 ) {    //S ingle RR rotation
                    tree.rre++;
                    p.right = p1.left;
                    p1.left = p;
                    if ( p1.bal == 0 ) {
                        p.bal = 1;
                        p1.bal = -1;
                        tree.changed = false;
                    } else {
                        p.bal = 0;
                        p1.bal = 0;
                    }
                    p = p1;
                } else {				// Double RL rotation
                    tree.rle++;
                    AvlNode p2 = p1.left;
                    p1.left = p2.right;
                    p2.right = p1;
                    p.right = p2.left;
                    p2.left = p;
                    if ( p2.bal == 1 ) {
                        p.bal = -1;
                    } else {
                        p.bal = 0;
                    }
                    if ( p2.bal == -1 ) {
                        p1.bal = 1;
                    } else {
                        p1.bal = 0;
                    }
                    p = p2;
                    p.bal = 0;
                }
                break;
            default:
                {
                    throw new RuntimeException("\n" + "bal = " + p.bal + " is out of range in AvlNode.balanceLeft\n");
                }
        }
        return p; // The root of the rebalanced sub-tree
    }
    
    /**
     * <p>
     * The {@code balanceRight} method rebalances following deletion of a right AvlNode.
     * 
     * This method is adapted from the {@code balance2} procedure on pp. 224-225 of Nickaus
     * Wirth's textbook, "Algorithms + Data Structures = Programs" (Prentice-Hall, 1976).
     * 
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    private AvlNode balanceRight(final AvlTree tree) {
        
        AvlNode p = this;
        
        switch ( p.bal ) {
            case 1:                     // Balance restored
                p.bal = 0;
                break;
            case 0:                     // Tree has become more unbalanced
                p.bal = -1;
                tree.changed = false;
                break;
            case -1:                    // Tree must be rebalanced
                AvlNode p1 = p.left;
                if ( p1.bal <= 0 ) {    // Single LL rotation
                    tree.lle++;
                    p.left = p1.right;
                    p1.right = p;
                    if ( p1.bal == 0 ) {
                        p.bal = -1;
                        p1.bal = 1;
                        tree.changed = false;
                    } else {
                        p.bal = 0;
                        p1.bal = 0;
                    }
                    p = p1;
                } else {				// Double LR rotation
                    tree.lre++;
                    AvlNode p2 = p1.right;
                    p1.right = p2.left;
                    p2.left = p1;
                    p.left = p2.right;
                    p2.right = p;
                    if ( p2.bal == -1 ) {
                        p.bal = 1;
                    } else {
                        p.bal = 0;
                    }
                    if ( p2.bal == 1 ) {
                        p1.bal = -1;
                    } else {
                        p1.bal = 0;
                    }
                    p = p2;
                    p.bal = 0;
                }
                break;
            default:
                {
                    throw new RuntimeException("\n" + "bal = " + p.bal + " is out of range in AvlNode.balanceRight\n");
                }
        }
        return p; // The root of the rebalanced sub-tree
    }
    
    /**
     * <p>
     * The {@code eraseLeft} method over-writes the fields of the AvlNode
     * to be deleted with the fields of the leftmost AvlNode of the right
     * sub-tree, and then removes that leftmost node from the sub-tree.
     * 
     * This method is adapted from the {@code del} procedure on p. 225 of
     * Nicklaus Wirth's textbook, "Algorithms + Data Structures = Programs"
     * (Prentice-Hall, 1976).
     * 
     * @param q - the {@code AvlNode} whose fields are to be replaced
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    AvlNode eraseLeft(AvlNode q,
                      final AvlTree tree)
    {
        
        AvlNode p = this;
        
        if ( p.left != null ) {
            p.left = left.eraseLeft( q, tree );
            if ( tree.changed == true ) {
                p = balanceLeft( tree );   // Height has changed, so rebalance sub-tree
            }
        } else {
            tree.removedNode = p;          // AvlNode p will be removed from the tree...
            final long[] tmpT = q.tuple;   // ...so swap the fields of p and q
            final KdNode tmpK = q.kdTreeNode;
            final short tmpI = q.kdTreeIndex;
            q.tuple = p.tuple;
            q.kdTreeNode = p.kdTreeNode;
            q.kdTreeIndex = p.kdTreeIndex;
            if (q.kdTreeNode != null) {
                q.kdTreeNode.avlTreeNode = q;
            }
            p.tuple = tmpT;
            p.kdTreeNode = tmpK;           // tmpK.values is empty for subsequent KdNode.erase
            p.kdTreeIndex = tmpI;
            if (p.kdTreeNode != null) {
                p.kdTreeNode.avlTreeNode = p;
            }
            p = p.right;                   // Replace AvlNode p with right branch...
            tree.changed = true;           // ...and signal that height has changed
        }
        return p; // The root of the sub-tree
    }
    
    /**
     * <p>
     * The {@code eraseRight} method over-writes the fields of the AvlNode
     * to be deleted with the fields of the rightmost AvlNode of the left
     * sub-tree, and then removes that rightmost node from the sub-tree.
     * 
     * This method is adapted from the {@code del} procedure on p. 225 of
     * Nicklaus Wirth's textbook, "Algorithms + Data Structures = Programs"
     * (Prentice-Hall, 1976).
     * 
     * @param q - the {@code AvlNode} whose fields are to be replaced
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    private AvlNode eraseRight(AvlNode q,
                               final AvlTree tree)
    {
        
        AvlNode p = this;
        
        if ( p.right != null ) {
            p.right = p.right.eraseRight( q, tree );
            if ( tree.changed == true ) {
                p = balanceRight( tree );  // Height has changed, so rebalance sub-tree
            }
        } else {
            tree.removedNode = p;          // AvlNode p will be removed from the tree...
            final long[] tmpT = q.tuple;   // ...so swap the fields of p and q
            final KdNode tmpK = q.kdTreeNode;
            final short tmpI = q.kdTreeIndex;
            q.tuple = p.tuple;
            q.kdTreeNode = p.kdTreeNode;
            q.kdTreeIndex = p.kdTreeIndex;
            if (q.kdTreeNode != null) {
                q.kdTreeNode.avlTreeNode = q;
            }
            p.tuple = tmpT;
            p.kdTreeNode = tmpK;           // tmpK.values is empty for subsequent KdNode.erase
            p.kdTreeIndex = tmpI;
            if (p.kdTreeNode != null) {
                p.kdTreeNode.avlTreeNode = p;
            }
            p = p.left;                    // Replace AvlNode p with left branch...
            tree.changed = true;           // ...and signal that height has changed
        }
        return p; // T the root of the sub-tree
    }
    
    /**
     * <p>
     * The {@code erase} method searches the AVL tree for the existence of
     * of a super-key. If the super-key is found, and the {@code TreeSet} of
     * the AvlNode that contains the super-key contains the value, the value
     * is removed from the {@code TreeSet}. If the {@code TreeSet} is empty
     * following removal of the value, the AvlNode is removed from the tree,
     * AvlTree.removedNode is assigned a reference to the removed AvlNode,
     * and AvlTree.nodeCount is decremented.
     * 
     * This method is adapted from the {@code delete} procedure on p. 225 of
     * Nicklaus Wirth's textbook, "Algorithms + Data Structures = Programs"
     * (Prentice-Hall, 1976).
     * 
     * NOTE: the super-key does not permute cyclically because an AVL tree is unidimensional.
     * 
     * @param key - a {@code long}[]
     * @param value - a {@code String}
     * @param tree - an {@code AvlTree}
     * @return the root of the re-balanced sub-tree
     * </p>
     */
    protected AvlNode erase(final long[] key,
                            final String value,
                            final AvlTree tree)
    {
        AvlNode p = this;
        
        if ( MergeSort.superKeyCompare(key, p.tuple, 0) < 0 ) { // Search left branch?
            if ( p.left != null ) {
                p.left = p.left.erase( key, value, tree );
                if ( tree.changed ) {
                    p = balanceLeft( tree );     // Height has changed, so rebalance sub-tree
                }
            } else {
                tree.erased = false;             // Key is not in the tree, so no value is,,,
                tree.changed = false;            // ...erased and height has not changed
            }
        } else if ( MergeSort.superKeyCompare(key, p.tuple, 0) > 0 ) { // Search right branch?
            if ( p.right != null ) {
                p.right = p.right.erase( key, value, tree );
                if ( tree.changed ) {
                    p = balanceRight( tree );    // Height has changed, so rebalance sub-tree
                }
            } else {
                tree.erased = false;             // Key is not in the tree, so no value is...
                tree.changed = false;            // ...erased and height has not changed
            }
        } else {                                // Key is in the tree at this node...
            if (p.kdTreeNode != null) {         // so if kdTreeNode has been assigned...
                p.kdTreeNode.values.remove(value);  // ...erase the value from kdTreeNode.values...
            }
            tree.erased = true;                 // ...and signal erasure of that value
            if ( p.kdTreeNode == null || p.kdTreeNode.values.isEmpty() ) {
                // Either kdTreeNode hasn't been assigned, or kdTreeNode.values is empty, so either
                // replace AvlNode p with its child, or swap its fields with its predecessor or successor
                if ( p.right == null ) {
                    tree.removedNode = p;       // Right branch is null..
                    p = p.left;                 // ...so replace AvlNode p with its left child...
                    tree.changed = true;        // ...and signal that the height has changed
                } else if ( p.left == null ) {  
                    tree.removedNode = p;       // Left branch is null...
                    p = p.right;                // ...so replace AvlNode p with its right child...
                    tree.changed = true;        // ...and signal that the height has changed
                } else {
                    switch ( p.bal ) {          // Neither branch is null...
                        case 0: case -1:        // ...and left or neither branch is deeper, so swap...
                            p.left = p.left.eraseRight( p, tree ); // ...AvlNode p's fields with predecessor
                            if ( tree.changed == true ) {
                                p = balanceLeft( tree );
                            }
                            break;
                        case 1:                 // Right branch is deeper, so swap...
                            p.right = p.right.eraseLeft( p, tree ); // ...AvlNode p's fields with successor
                            if ( tree.changed == true ) {
                                p = balanceRight( tree );
                            }
                            break;
                        default:
                            {
                                throw new RuntimeException("\n" + "bal = " + p.bal +
                                                           " is out of range in AvlNode.erase\n");
                            }
                    }
                }
            }
        }
        return p; // The root of the rebalanced sub-tree
    }

    /**
     * <p>
     * The {@code verifyAvlTree} method computes the maximum height at each (@code AvlNode}
     * in the AvlTree and then checks the balance and sorted order at those {@code AvlNodes}.
     * 
     * @return the computed height
     * </p>
     */
    protected int verifyAvlTree()
    {
        // Compute the maximum height at each AvlNode.
        computeHeights();

        // Verify correct height and sorted order at each AvlNode.
        return verifyAvlNodes();
    }

    /**
     * <p>
     * The {@code verifyAvlNodes} method counts each {@code AvlNode},
     * and checks the balance and sorted order at that {@code AvlNode}.
     * 
     * NOTE: the super-key does not permute cyclically because an AVL tree is unidimensional.
     * 
     * @return the number of nodes in the AvlTree
     * </p>
     */
    private int verifyAvlNodes()
    {
        // Count this node and check its balance.
        int count = 1;
        if (bal != right.height - left.height) {
            throw new RuntimeException("\nbal = " + bal + "  right.height = " + right.height +
                                       "  left.height = " + left.height + "\n");
        }

        if (left != null) {
            if ( MergeSort.superKeyCompare(left.tuple, tuple, 0) >= 0) {
                throw new RuntimeException("\nleft >= this\n");
            }
            count += left.verifyAvlNodes();
        }

        if (right != null) {
            if ( MergeSort.superKeyCompare(right.tuple, tuple, 0) <= 0) {
                throw new RuntimeException("\nright <= this\n");
            }
            count += right.verifyAvlNodes();
        }

        return count;
    }

     /**
     * <p>
     * The {@code computeHeights} method computes the maximum height at each (@code AvlNode}.
     * </p>
     */
   private void computeHeights()
    {
        if (left != null) {
            left.computeHeights();
        }
        if (right != null) {
            right.computeHeights();
        }
        height = computeHeight();
    }
    
    /**
     * <p>
     * The {@code computeHeight} method computes the maximum height at an (@code AvlNode}.
     * 
     * @return the maximum height
     * </p>
     */
    private int computeHeight()
    {
        return 1 +  ( (getHeight(left) > getHeight(right) )
                     ? getHeight(left) : getHeight(right) ); 
    }

    /**
     * <p
     * The {@code getHeight} method returns the sub-tree height at an (@code AvlNode}.
     * 
     * @return the height
     * </p>
     */
    private int getHeight(final AvlNode node)
    {
        if (node == null) {
            return 0;
        } else {
            return node.height;
        }
    }
    
} // class AvlNode
