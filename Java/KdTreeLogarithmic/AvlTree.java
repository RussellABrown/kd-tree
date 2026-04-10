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
 * The {@code AvlNode} class stores a point of any number of dimensions
 * as well as references to the "less than" and "greater than" sub-trees,
 * a reference to a KdNode instance, a duplicates count, a balance, and
 * an index into the KdTreeLogarithmic.kdTrees array.
 * </p>
 */
public class AvlTree {
    
    private long nodeCount;
    protected boolean inserted, erased, changed; // "changed" means height changed
    protected long lle, lre, rle, rre, lli, lri, rli, rri; // rotation counters
    protected AvlNode root, insertedNode, removedNode;

    /**
     * <p>
     * AvlTree constructor
     * <p>
     */
    public AvlTree()
    {
        inserted = erased = changed = false;
        lle = lre = rle = rre = lli = lri = rli = rri = nodeCount = 0L;
        root = insertedNode = removedNode = null;
    }

   /**
     * <p>
     * The {@code contains} method searches the AVL tree for the existence
     * of a (key, value) pair, and returns true if the AVL tree contains it.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code true} if the tree contains the key; otherwise, {@code false}
     * </p>
     */
    protected boolean contains(final Pair coordinate)
    {
        if (root != null) {
            return root.contains( coordinate.getKey(), coordinate.getValue() );
        } else {
            return false;
        }
    }

    /**
     * <p>
     * The {@code find} method searches the AVL tree for the existence
     * of a (key, value) pair, and returns the AvlNode instance that contains it.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code AvlNode} reference if the tree contains the key; otherwise, {@code null}
     * </p>
     */
    protected AvlNode find(final Pair coordinate)
    {
        if (root != null) {
            return root.find( coordinate.getKey(), coordinate.getValue() );
        } else {
            return null;
        }
    }

    /**
     * <p>
     * The {@code insert} method searches the AVL tree for the existence
     * of a super-key, inserts the super-key and its associated value
     * and increments AvlTree.nodeCount if the super-key was not found,
     * but inserts only the associated value if the super-key was found.
     * 
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code true} if an AvlNode is inserted; otherwise {@code false}
     * </p>
     */
    protected boolean insert(final Pair coordinate)
    {
        return insert( coordinate.getKey(), coordinate.getValue() );
    }

    /**
     * <p>
     * The {@code insert} method searches the AVL tree for the existence
     * of a super-key, inserts the super-key and its associated value
     * and increments AvlTree.nodeCount if the super-key was not found,
     * but inserts only the associated value if the super-key was found.
     * 
     * @param key - a {@code long}[]
     * @param value - a {@code String}
     * @return {@code true} if an AvlNode is inserted; otherwise {@code false}
     * </p>
     */
    protected boolean insert(final long[] key,
                             final String value)
    {
        inserted = changed = false;
        insertedNode = null;
        if (root != null) {
            root = root.insert(key, value, this);
            if (insertedNode != null) {
                ++nodeCount;
            }
        } else {
            root = insertedNode = new AvlNode(key);
            inserted = true;
            ++nodeCount;
        }
        return inserted;
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
     * @param coordinate - a {@code Pair}<{@code long}[],{@code String}>
     * @return {@code true} if a value is removed; otherwise {@code false}
     * </p>
     */
    protected boolean erase(final Pair coordinate)
    {
        return erase( coordinate.getKey(), coordinate.getValue() );
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
     * @param key - a {@code long}[]
     * @param value - a {@code String}
     * @return {@code true} if a value is removed; otherwise {@code false}
     * </p>
     */
    protected boolean erase(final long[] key,
                            final String value)
    {
        erased = changed = false;
        removedNode = null;
        if (root != null) {
            root = root.erase(key, value, this);
            if (removedNode != null) {
                --nodeCount;
            }
        }
        return erased;
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
        if (root != null) {
            return root.verifyAvlTree();
        } else {
            return 0;
        }
    }

    /**
     * <p>
     * The {@code size} method returns the size of the {@code AvlTree}.
     * 
     * @return the tree size
     * </p>
     */
    protected long size()
    {
        return nodeCount;
    }

    /**
     * <p>
     * The {@code size} method returns the size of the {@code AvlTree}.
     * 
     * @return the tree size
     * </p>
     */
    protected boolean isEmpty()
    {
        return (size() == 0);
    }

    /**
     * <p>
     * The {@code getHeight} method returns the maximum height of the {@code AvlTree}.
     * 
     * @return the height at the root {@code AvlNode}
     * </p>
     */
    protected int getHeight()
    {
        if (root != null) {
            return root.height;
        } else {
            return 0;
        }
    }

} // class AvlTree
