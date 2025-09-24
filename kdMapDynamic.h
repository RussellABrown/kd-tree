/*
 * Modifications Copyright (c) 2025 Russell A. Brown
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

#ifndef KD_MAP_DYNAMIC_H
#define KD_MAP_DYNAMIC_H

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::cout;
using std::endl;
using std::ostringstream;
using std::pair;
using std::runtime_error;
using std::set;
using std::streamsize;
using std::vector;

/* The maximum allowed height difference between a node's < and > subtrees. */
#ifndef HEIGHT_DIFF
#define HEIGHT_DIFF (1)
#endif

/* The size of the histogram vectors. */
#ifndef HISTOGRAM_SIZE
#define HISTOGRAM_SIZE (25)
#endif

/* The size of the maximum value vectors. */
#ifndef MAXIMUM_SIZE
#define MAXIMUM_SIZE (25)
#endif

/* Convert microseconds to seconds for use with std::chrono */
static const double MICROSECONDS_TO_SECONDS = 1000000.;

/* This type is the signed equivalent of size_t and might be equivalent to intmax_t */
typedef streamsize signed_size_t;

/*
 * Forward references to all classes to support any order of compilation
 *
 * However, this kdMapDynamic.h file must be compiled before kdMapNode.h,
 * kdMapKnlogn.h, and kdMapNlogn.h because those files contain #if and #ifdef
 * directives that check whether KD_MAP_DYNAMIC_H is defined.
 */
template <typename, typename>
class KdTreeDynamic;

template <typename, typename>
class KdTree;

template <typename, typename>
class KdNode;

template <typename, typename>
class MergeSort;

template <typename, typename>
class NearestNeighborHeap;

/*
 * The KdTreeDynamic class is derived from the KdTree class and
 * accessess KdTree::root so that an instance of KdTreeDynamic
 * may use KdTree member functions that require access to root.
 */
template <typename K, typename V>
class KdTreeDynamic : public KdTree<K,V>
{
private:
    signed_size_t maxSubmitDepth = -1;
    size_t cutoff;
    bool inserted = false, erased = false, changed = false;

    /*
     * This is the basic constructor.
     *
     * Calling parameters:
     * 
     * maxSubmitDepth (IN) the maximum tree depth for creating a child thread
     * cutoff (IN) the minimum size of a subtree to rebuild via multiple threads
     */
public:
    KdTreeDynamic(signed_size_t const maxSubmitDepth,
                  size_t const cutoff) {

        this->maxSubmitDepth = maxSubmitDepth;
        this->cutoff = cutoff;
    }

    /*
     * This constructor assigns the KdTree::root node.
     *
     * Calling parameters:
     * 
     * maxSubmitDepth (IN) the maximum tree depth for creating a child thread
     * cutoff (IN) the minimum size of a subtree to rebuild via multiple threads
     * root (IN) the KdTree::root node
     */
public:
    KdTreeDynamic(signed_size_t const maxSubmitDepth,
                  size_t const cutoff,
                  KdNode<K,V>* const root) {

        this->maxSubmitDepth = maxSubmitDepth;
        this->cutoff = cutoff;
        KdTree<K,V>::root = root;
    }

    /*
     * If KD_MAP_DYNAMIC_H is defined, the ~KdTree destructor does not
     * delete KdTree::root to prevent recursive deletion of the k-d tree
     * whose persistence is required by KdTreeDynamic::rebuildSubTree, so
     * so delete KdTree::root here.
     */
public:
    ~KdTreeDynamic() {
        delete KdTree<K,V>::root;
    }

    /*
     * Search the tree for the existence of a key.
     *
     * Calling parameter:
     *
     * @param coordinate (IN) the key-value pair to search for
     * 
     * @return true if the key was found; otherwise, false
     */
public:
    inline bool contains(pair<vector<K>, V> const& coordinate) {

        if (KdTree<K,V>::root != nullptr) {
            vector<K> const tuple = coordinate.first;
            signed_size_t const dim = tuple.size();
            K* const key = const_cast<K* const>(tuple.data());
            V const value = coordinate.second;
            return contains(KdTree<K,V>::root, key, value, dim, 0);
        } else {
            return false;
        }
    }
    
   /* Search the tree for the existence of a key.
    *
    * Calling parameters:
    *
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple to search for
     * @param value (IN) the value to search for
     * @param dim (IN) the number of dimensions
     * @param partition (IN) the leading dimension that is permuted cyclically
    * 
    * @return true if the key was found; otherwise, false
    */
private:
    inline bool contains(KdNode<K,V>* const node,
                         K* const key,
                         V const& value,
                         signed_size_t const dim,
                         signed_size_t p) {

        KdNode<K,V>* ptr = node;
        while ( ptr != nullptr ) {
            if (MergeSort<K,V>::superKeyCompare(key, ptr->tuple, p, dim) < 0) {
                ptr = ptr->ltChild;
            } else if (MergeSort<K,V>::superKeyCompare(key, ptr->tuple, p, dim) > 0) {
                ptr = ptr->gtChild;
            } else {
                // found the key, so check for the value
                return ptr->values->contains(value);
            }
            // Permute the most significant dimension p cyclically using
            // a fast alternative to the modulus operator for p <= dim.
            ++p;
            p = (p < dim) ? p : 0;
        }
        return false; // didn't find the key
    }
    
    /*
     * Insert a key into the tree.
     *
     * Calling parameter:
     *
     * @param coordinate (IN) the key-value pair to insert
     * 
     * This function modifies KdTreeDynamic::root,
     * KdTreeDynamic::inserted, and KdTreeDynamic::changed.
     * 
     * @return true if the key was inserted into a node's set; otherwise, false
     */
public:
    inline bool insert(pair<vector<K>, V> const& coordinate) {

        inserted = changed = false;
        vector<K> const tuple = coordinate.first;
        signed_size_t const dim = tuple.size();
        K* const key = const_cast<K* const>(tuple.data());
        V const value = coordinate.second;
        if (KdTree<K,V>::root != nullptr) {
            KdTree<K,V>::root = insert(KdTree<K,V>::root, key, value, dim, 0);

            // Has the height changed due to an insertion?
            if (changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(KdTree<K,V>::root) ) {
                    KdTree<K,V>::root->height = computeHeight(KdTree<K,V>::root);
                } else {
                    KdTree<K,V>::root = rebuildSubTree(KdTree<K,V>::root, dim, 0);
                }
            }
        } else {
            KdTree<K,V>::root = new KdNode<K,V>();
            KdTree<K,V>::root->height = 1;
            // The tuple member field will be deleted by the ~KdNode destructor.
            KdTree<K,V>::root->tuple = new K[dim];
            for (signed_size_t i = 0; i < dim; ++i) {
                KdTree<K,V>::root->tuple[i] = key[i];
            }
            // The values member field will be cleared by the ~KdNode destructor.
            KdTree<K,V>::root->values->insert(value);
            inserted = changed = true;
        }
        return inserted;
    }

    /*
     * Insert a key into the subtree.
     *
     * Calling parameters:
     *
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple to insert
     * @param value (IN) the value to insert
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension that is permuted cyclically
     * 
     * @return the root of the possibly rebalanced subtree
     */
private:
    KdNode<K,V>* insert(KdNode<K,V>* const node,
                        K* const key,
                        V const& value,
                        signed_size_t dim,
                        signed_size_t p) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // The return value because the node argument to this function is read only.
        KdNode<K,V>* nodePtr = node;

        // Determine which child to search recursively for an insertion point.
        if (MergeSort<K,V>::superKeyCompare(key, nodePtr->tuple, p, dim) < 0) {
            if (nodePtr->ltChild != nullptr) {
                nodePtr->ltChild = insert(nodePtr->ltChild, key, value, dim, p+1);
            } else {
                // The tree does not contain the key, so insert it.
                nodePtr->ltChild = new KdNode<K,V>();
                nodePtr->ltChild->height = 1;
                // The tuple member field will be deleted by the ~KdNode destructor.
                nodePtr->ltChild->tuple = new K[dim];
                for (signed_size_t i = 0; i < dim; ++i) {
                    nodePtr->ltChild->tuple[i] = key[i];
                }
                // The values member field will be cleared by the ~KdNode destructor.
                nodePtr->ltChild->values->insert(value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else if (MergeSort<K,V>::superKeyCompare(key, nodePtr->tuple, p, dim) > 0) {
            if (nodePtr->gtChild != nullptr) {
                nodePtr->gtChild = insert(nodePtr->gtChild, key, value, dim, p+1);
            } else {
                // The tree does not contain the key, so insert it.
                nodePtr->gtChild = new KdNode<K,V>();
                nodePtr->gtChild->height = 1;
                // The tuple member field will be deleted by the ~KdNode destructor.
                nodePtr->gtChild->tuple = new K[dim];
                for (signed_size_t i = 0; i < dim; ++i) {
                    nodePtr->gtChild->tuple[i] = key[i];
                }
                // The values member field will be cleared by the ~KdNode destructor.
                nodePtr->gtChild->values->insert(value);
                // The value was inserted and the tree height changed.
                inserted = changed = true;
            }
        } else {
            // The tree contains the key, so insert the value,
            // even if the values set might already contain
            // it. The tree height doesn't change.
            nodePtr->values->insert(value);
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
                nodePtr->height = computeHeight(nodePtr);
            } else {
                nodePtr = rebuildSubTree(nodePtr, dim, p);
            }
       }
        return nodePtr;
    }

    /*
     * Erase a key from the tree.
     *
     * Calling parameter:
     * 
     * @param coordinate (IN) the key-value pair to erase
     * 
     * @return true if the key erased from a node's set; otherwise, false
     */
public:
    inline bool erase(pair<vector<K>, V> const& coordinate) {

        erased = changed = false;
        if (KdTree<K,V>::root != nullptr) {
            vector<K> const tuple = coordinate.first;
            signed_size_t const dim = tuple.size();
            K* const key = const_cast<K* const>(tuple.data());
            V const value = coordinate.second;
            KdTree<K,V>::root = erase(KdTree<K,V>::root, key, value, dim, 0);

            // Has the height changed due to an erasure?
            if (KdTree<K,V>::root != nullptr && changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(KdTree<K,V>::root) ) {
                    KdTree<K,V>::root->height = computeHeight(KdTree<K,V>::root);
                } else {
                    KdTree<K,V>::root = rebuildSubTree(KdTree<K,V>::root, dim, 0);
                }
            }
        }
        return erased;
    }

    /*
     * Erase a node from the subtree.
     * 
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple to erase
     * @param value (IN) the value to erase
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension that is permuted cyclically
     * @param clearSet (IN) clear all values in the set (default false)
     * 
     * @return the root of the possibly rebalanced subtree
     */
private:
    KdNode<K,V>* erase(KdNode<K,V>* const node,
                      K* const key,
                      V const& value,
                      signed_size_t const dim,
                      signed_size_t p,
                      bool const clearSet = false) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // The return value because the node argument to this function is read only.
        KdNode<K,V>* nodePtr = node;

        // Determine which child to search recursively for a deletion point.
        if (MergeSort<K,V>::superKeyCompare(key, nodePtr->tuple, p, dim) < 0) {
            if (nodePtr->ltChild != nullptr) {
                nodePtr->ltChild = erase(nodePtr->ltChild, key, value, dim, p+1, clearSet);

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
                        nodePtr->height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, dim, p);
                    }
                }
            } else {
                // The tree does not contain the key.
                erased = changed = false;
            }
        } else if (MergeSort<K,V>::superKeyCompare(key, nodePtr->tuple, p, dim) > 0) {
            if (nodePtr->gtChild != nullptr) {
                nodePtr->gtChild = erase(nodePtr->gtChild, key, value, dim, p+1, clearSet);
             
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
                        nodePtr->height = computeHeight(nodePtr);
                    } else {
                        nodePtr = rebuildSubTree(nodePtr, dim, p);
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
            if ( !clearSet && !nodePtr->values->contains(value) ) {
                erased = changed = false;
            } else {
                // Either clearSet is true, or the values set
                // contains the value. So if clearSet is true,
                // clear the values set; otherwise, remove the
                // value from the values set.
                if (clearSet) {
                    nodePtr->values->clear();
                } else {
                    nodePtr->values->erase(value);
                }

                // Is the values set now empty?
                if ( !nodePtr->values->empty() ) {
                    // The values set is not empty, so a value was erased
                    // but the tree height has not changed.
                    erased = true;
                    changed = false;
                } else {
                    // The values set is empty, so a value was erased and
                    // the tree height will change via removal of a node.
                    erased = changed = true;

                    // Does the node have only a < child?
                    if (nodePtr->ltChild != nullptr && nodePtr->gtChild == nullptr) {
                        // Yes, the node has only a < child, so does the subtree
                        // rooted at the < child contain <= 3 nodes?
                        //
                        // Checking the height prior to checking the countNodes result
                        // avoids a time-consuming call of countNodes for a large subtree.

#ifdef ENABLE_1TO3

                        size_t nodeCount;
                        if (nodePtr->ltChild->height <= 3
                            && (nodeCount = countNodes(nodePtr->ltChild)) <= 3)
                        {
                            // Yes, the < child subtree contains <= 3 nodes. So, rebuild
                            // the < child subtree using the leading dimension p at this
                            // level of the tree, and delete the node. This approach
                            // avoids the need to find a predecessor node.
                            //
                            // Note that rebuildSubTree1to3 computes the height
                            // of the rebuilt < child subtree.
                            KdNode<K,V>* const tempPtr = nodePtr;
                            nodePtr = rebuildSubTree1to3(nodePtr->ltChild, nodeCount, dim, p);
                            tempPtr->ltChild = nullptr; // Prevent recursive deletion.
                            delete tempPtr;
                        } else

#endif // ENABLE_1TO3

                        {
                            // No, the < child subtree contains > 3 nodes. So, find
                            // the immediate predecessor node, copy the tuple and the
                            // values set from that predecessor node to the one-child
                            // node, delete the predecessor node recursively (clearing
                            // its values set), and recompute the height along the path
                            // back to the < child, including that child.
                            KdNode<K,V>* predecessor = nodePtr->ltChild;
                            predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                            for (signed_size_t i = 0; i < dim; ++i) {
                                nodePtr->tuple[i] = predecessor->tuple[i];
                            }
                            nodePtr->values->insert(predecessor->values->begin(), predecessor->values->end());
                            // value is a dummy argument because the clearSet argument is true
                            nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, value, dim, p+1, true);

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
                                nodePtr->height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, dim, p);
                            }
                        }
                    }
                    // Does the node have only a > child?
                    else if (nodePtr->gtChild != nullptr && nodePtr->ltChild == nullptr) {
                        // Yes, the node has only a > child, so does the subtree
                        // rooted at the > child contain <= 3 nodes?
                        //
                        // Checking the height prior to checking the countNodes result
                        // avoids a time-consuming call of countNodes for a large subtree.

#ifdef ENABLE_1TO3

                        size_t nodeCount;
                        if (nodePtr->gtChild->height <= 3
                            && (nodeCount = countNodes(nodePtr->gtChild)) <= 3)
                        {
                            // Yes, the > child subtree contains <= 3 nodes. So, rebuild
                            // the > child subtree using the leading dimension p at this
                            // level of the tree, and delete the node. This approach
                            // avoids the need to find a predecessor node.
                            //
                            // Note that rebuildSubTree1to3 computes the height
                            // of the rebuilt > child subtree.
                            KdNode<K,V>* const tempPtr = nodePtr;
                            nodePtr = rebuildSubTree1to3(nodePtr->gtChild, nodeCount, dim, p);
                            tempPtr->gtChild = nullptr; // Prevent recursive deletion.
                            delete tempPtr;
                        } else

#endif // ENABLE_1TO3

                        {
                            // No, the > child subtree contains > 3 nodes. So, find
                            // the immediate successor node, copy the tuple and the
                            // values set from that successor node to the one-child
                            // node, delete the successor node recursively (clearing
                            // its values set), and recompute the height along the path
                            // back to the > child, including that child.
                           KdNode<K,V>* successor = nodePtr->gtChild;
                            successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                            for (signed_size_t i = 0; i < dim; ++i) {
                                nodePtr->tuple[i] = successor->tuple[i];
                            }
                            nodePtr->values->insert(successor->values->begin(), successor->values->end());
                            // value is a dummy argument because the clearSet argument is true
                            nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, value, dim, p+1, true);

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
                                nodePtr->height = computeHeight(nodePtr);
                            } else {
                                nodePtr = rebuildSubTree(nodePtr, dim, p);
                            }
                        }
                    }
                    // If the node has no children, delete the node.
                    else if (nodePtr->ltChild == nullptr && nodePtr->gtChild == nullptr) {
                        KdNode<K,V>* const tempPtr = nodePtr;
                        nodePtr = nullptr;
                        delete tempPtr;
                    }
                    // The node has two children.
                    else { 
                        // So does the subtree rooted at this two-child node 
                        // have <= 3 nodes, excluding this two-child node?
                        //
                        // Checking the height prior to checking the countNodesSkipRoot result
                        // avoids a time-consuming call of countNodesSkipRoot for a large subtree.

#ifdef ENABLE_1TO3

                        size_t nodeCount;
                        if (nodePtr->height <= 3
                            && (nodeCount = countNodesSkipRoot(nodePtr)) <= 3)
                        {
                            // Yes, the subtree rooted at this two-child node contains
                            // <= 3 nodes, excluding this two-child node, so rebuild
                            // the subtree excluding this two-child node and delete
                            // this two-child node. This approach avoids the need to
                            // find a predecessor or successor node.
                            KdNode<K,V>* const tempPtr = nodePtr;
                            nodePtr = rebuildSubTreeSkipRoot1to3(nodePtr, nodeCount, dim, p);
                            tempPtr->ltChild = tempPtr->gtChild = nullptr; // Prevent recursive deletion.
                            delete tempPtr;
                        } else

#endif // ENABLE_1TO3

                        {
                            // No, the subtree rooted at this two-child node contains
                            // > 3 nodes, excluding this two-child node, so replace this
                            // two-child node by either its predecessor or its successor.
                            //
                            // If ENABLE_PREFERRED_TEST is defined, select the replacement
                            // node from the taller of the child subtrees.

#ifdef ENABLE_PREFERRED_TEST

                            if ( getHeight(nodePtr->ltChild) >= getHeight(nodePtr->gtChild) )
                            {
                                // Find the node with the largest super-key in the
                                // subtree rooted at the < child, which is the
                                // predecessor node. Copy the predecessor node's tuple
                                // and values set to this two-child node, delete the
                                // predecessor node recursively (clearing its values set),
                                // and recompute the heights along the path from the
                                // predecessor node to (but excluding) this two-child node.
                                KdNode<K,V>* predecessor = nodePtr->ltChild;
                                predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                                for (signed_size_t i = 0; i < dim; ++i) {
                                    nodePtr->tuple[i] = predecessor->tuple[i];
                                }
                                nodePtr->values->insert(predecessor->values->begin(), predecessor->values->end());
                                // value is a dummy argument because the clearSet argument is true
                                nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, value, dim, p+1, true);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( isBalanced(nodePtr) ) {
                                    nodePtr->height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, dim, p);
                                }
                            } else

#endif // ENABLE_PREFERRED_TEST

                            {
                                // Find the node with the smallest super-key in the
                                // subtree rooted at the > child, which is the
                                // successor node. Copy the successor node's tuple
                                // and values set to this two-child node, delete the
                                // successor node recursively (clearing its values set),
                                // and recompute the heights along the path from the
                                // successor node to (but excluding) this two-child node.
                                KdNode<K,V>* successor = nodePtr->gtChild;
                                successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                                for (signed_size_t i = 0; i < dim; ++i) {
                                    nodePtr->tuple[i] = successor->tuple[i];
                                }
                                nodePtr->values->insert(successor->values->begin(), successor->values->end());
                                // value is a dummy argument because the clearSet argument is true
                                nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, value, dim, p+1, true);

                                // The height may have changed, so if the subtree
                                // rooted at this two-child node remains balanced,
                                // compute the height at this node; otherwise,
                                // rebuild the subtree to rebalance it, which also
                                // computes the height. Because rebuilding the subtree
                                // recycles its nodes, the node argument to this
                                // erase function might no longer specify the root
                                // of the subtree.
                                if ( isBalanced(nodePtr) ) {
                                    nodePtr->height = computeHeight(nodePtr);
                                } else {
                                    nodePtr = rebuildSubTree(nodePtr, dim, p);
                                }
                            }
                        }
                    }
                }
            }
        }
        return nodePtr;
    }

    /*
     * Search the subtree for an immediate predecessor node,
     * which is the node with the largest super-key.
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param predecessor (MODIFIED) the current predecessor node
     * @param dim (IN) the number of dimensions
     * @param p0 (IN) the leading dimension of the node to be replaced
     * @param p (IN) the leading dimension that is permuted cyclically
     * 
     * @return the updated predecessor node
     */
private:
    KdNode<K,V>* findPredecessor(KdNode<K,V>* const node,
                                 KdNode<K,V>* const predecessor,
                                 signed_size_t const dim,
                                 signed_size_t const p0,
                                 signed_size_t p) {

        // The return value because the predecessor argument to this function is read only.
        KdNode<K,V>* pred = predecessor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // Does the leading dimension at this node equal
        // the leading dimension at the node to be replaced?
        if (p == p0) {
            // Yes, the leading dimensions are equal, so if this
            // node has a > child, follow that child; otherwise,
            // this node is a potential predecessor, so check it.
            if (node->gtChild != nullptr) {
                pred = findPredecessor(node->gtChild, pred, dim, p0, p+1);
            } else {
                pred = checkPredecessor(node, pred, dim, p0);
            }
        } else {
            // No, the leading dimensions are not equal; hence, this
            // node is a potential predecessor, so check it, and
            // then inspect both children.
            pred = checkPredecessor(node, pred, dim, p0);
            if (node->ltChild != nullptr) {
                pred = findPredecessor(node->ltChild, pred, dim, p0, p+1);
            }
            if (node->gtChild != nullptr) {
                pred = findPredecessor(node->gtChild, pred, dim, p0, p+1);
            }
        }
        return pred;
    }

    /*
     * Check a potential predecessor node to see whether it is
     * a larger predecessor than the current predecessor node.
     * 
     * Calling parameters:
     *
     * @param node (IN) a potential predecessor node
     * @param predecessor (MODIFIED) the current predecessor node
     * @param dim (IN) the number of dimensions
     * @param p0 (IN) the leading dimension of the node to be replaced
     * 
     * @return the updated predecessor node
     */
private:
    KdNode<K,V>* checkPredecessor(KdNode<K,V>* const node,
                                  KdNode<K,V>* const predecessor,
                                  signed_size_t const dim,
                                  signed_size_t const p0) {

        // The return value because the predecessor argument to this function is read only.
        KdNode<K,V>* pred = predecessor;

        // Update the current precedessor if the potential predecessor is larger.
        if (MergeSort<K,V>::superKeyCompare(node->tuple, pred->tuple, p0, dim) > 0) {
            pred = node;
        }
        return pred;
    }

    /*
     * Search the subtree for an immediate successor node,
     * which is the node with the smallest super-key.
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param successor (MODIFIED) the current successor node
     * @param dim (IN) the number of dimensions
     * @param p0 (IN) the leading dimension of the node to be replaced
     * @param p (IN) the leading dimension that is permuted cyclically
     * 
     * @return the updated successor node
     */
private:
    KdNode<K,V>* findSuccessor(KdNode<K,V>* const node,
                               KdNode<K,V>* const successor,
                               signed_size_t const dim,
                               signed_size_t const p0,
                               signed_size_t p) {

        // The return value because the successor argument to this function is read only.
        KdNode<K,V>* succ = successor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // Does the leading dimension at this node equal
        // the leading dimension at the node to be replaced?
         if (p == p0) {
            // Yes, the leading dimensions are equal, so if this
            // node has a < child, follow that child; otherwise,
            // this node is a potential successor, so check it.
            if (node->ltChild != nullptr) {
                succ = findSuccessor(node->ltChild, succ, dim, p0, p+1);
            } else {
                succ = checkSuccessor(node, succ,  dim, p0);
            }
        } else {
            // No, the leading dimensions are not equal; hence, this
            // node is a potential successor, so check it, and
            // then inspect both children
            succ = checkSuccessor(node, succ, dim, p0);
            if (node->ltChild != nullptr) {
                succ = findSuccessor(node->ltChild, succ, dim, p0, p+1);
            }
            if (node->gtChild != nullptr) {
                succ = findSuccessor(node->gtChild, succ, dim, p0, p+1);
            }
        }
        return succ;
    }

    /*
     * Check a potential successor node to see whether it is
     * a smaller successor than the current successor node.
     * 
     * Calling parameters:
     *
     * @param node (IN) a potential successor node
     * @param successor (MODIFIED) the current successor node
     * @param dim (IN) the number of dimensions
     * @param p0 (IN) the leading dimension of the node to be replaced
     * 
     * @return the updated successor node
     */
private:
    KdNode<K,V>* checkSuccessor(KdNode<K,V>* const node,
                                KdNode<K,V>* const successor,
                                signed_size_t const dim,
                                signed_size_t p0) {

        // The return value because the successor argument to this function is read only.
        KdNode<K,V>* succ = successor;

        // Update the current successor if the potential successor is smaller.
        if (MergeSort<K,V>::superKeyCompare(node->tuple, succ->tuple, p0, dim) < 0) {
            succ = node;
        }
        return succ;
    }

    /*
     * Rebuild a subtree to rebalance it.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the balanced subtree
     */
private:
    inline KdNode<K,V>* rebuildSubTree(KdNode<K,V>* const node,
                                       size_t const dim,
                                       signed_size_t const p) {

        // Walk the subtree to obtain a vector of pointers to KdNodes.
        // Count the nodes in the subtree, allocate a sized vector of KdNode
        // pointers, and assign to each vector element a pointer to each KdNode.
        size_t count = countNodes(node);
        vector<KdNode<K,V>*> kdNodes(count);
        size_t index = 0;
        getSubTree(node, kdNodes, index);

#ifndef ENABLE_1To3

        // Does the subtree contain 3 KdNodes or fewer?
        //
        // Note that there is no need for this test if
        // ENABLE_1TO3 is defined because, in that case,
        // rebuildSubTree1to3 is called for 3 nodes or fewer.
        if (count <= 3) {
            // Yes, so rebuild the subtree by explicitily comparing the
            // nodes' super-keys.
             return rebuildSubTree1to3(node, kdNodes, dim, p);
        } else

#endif
        
        {
            // No, the subtree contains more than 3 KdNodes,
            // so rebuild the subtree via KdTree::createKdTree.

            // Check the count and kdNodes.size() for consistency.
            if (count != kdNodes.size()) {
                ostringstream buffer;
                buffer << "count = " << count << "  size = " << kdNodes.size() << endl;
                throw runtime_error(buffer.str());
            }

            // Call KdTree::createKdTree to rebuild the subtree, which
            // invalidates the node argument to this rebuildSubTree function.
            //
            // Do not use mutiple threads to build the subtree, even
            // though multiple threads are available, unless the size of the
            // subtree is sufficiently large to justify spawning child threads.
            KdTree<K,V>* subTree;
            if (count < cutoff) {
                subTree =
                    KdTree<K,V>::createKdTree(kdNodes, dim, -1, p);
            } else {
                subTree =
                    KdTree<K,V>::createKdTree(kdNodes, dim, maxSubmitDepth, p);
            }

            // Delete the KdTree instance but not its root
            // (see the ~KdTree destructor) so that the 
            // rebuilt subtree will be added to the k-d tree.
            KdNode<K,V>* const subTreeRoot = subTree->root;
            delete subTree;
            return subTreeRoot;
        }
    }

    /*
     * Build a subtree that contains 1 to 3 KdNode instances.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param nodeCount (IN) the number of nodes in the subtree
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K,V>* rebuildSubTree1to3(KdNode<K,V>* const node,
                                           size_t const nodeCount,
                                           size_t const dim,
                                           signed_size_t const p) {

        // Allocate a sized vector of KdNode pointers, and walk the subtree
        // to assign to each vector element a pointer to a KdNode instance.
        vector<KdNode<K,V>*> kdNodes(nodeCount);
        size_t index = 0;
        getSubTree(node, kdNodes, index);
        return rebuildSubTree1to3(node, kdNodes, dim, p);
    }

    /*
     * Build a subtree that contains 1 to 3 KdNode instances.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param kdNodes (IN) a vector of kdNode pointers
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K,V>* rebuildSubTree1to3(KdNode<K,V>* const node,
                                           vector<KdNode<K,V>*> const& kdNodes,
                                           size_t const dim,
                                           signed_size_t const p) {

        // The return value because the node argument to this function is read only.
        KdNode<K,V>* ptr = node;

        if (kdNodes.size() == 1) {

            // The subtree contains 1 node, so it is the root of the rebuilt subtree.
            ptr = kdNodes[0];
            ptr->height = 1;

        } else if (kdNodes.size() == 2) {

            // The subtree contains 2 nodes, so the first node is the root of the rebuilt
            // subtree. Compare super-keys to determine whether second node is the < child
            // or the > child of the first node.
            ptr = kdNodes[0];
            ptr->height = 2;
            if (MergeSort<K,V>::superKeyCompare(kdNodes[0]->tuple, kdNodes[1]->tuple, p, dim) > 0) {
                ptr->ltChild = kdNodes[1];
                ptr->ltChild->height = 1;
            } else {
                ptr->gtChild = kdNodes[1];
                ptr->gtChild->height = 1;
            }

        } else if (kdNodes.size() == 3) {

            // The subtree contains 3 nodes, so compare their super-keys to determine which
            // node is the median node. The node with the smallest super-key is the < child.
            // The node with the larges super-key is the > child.
            if (MergeSort<K,V>::superKeyCompare(kdNodes[0]->tuple, kdNodes[1]->tuple, p, dim) < 0) {
                // kdNodes[0]->tuple < kdNodes[1]->tuple
                if (MergeSort<K,V>::superKeyCompare(kdNodes[1]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                    // kdNodes[0]->tuple < kdNodes[1]->tuple < kdNodes[2]->tuple
                    ptr = kdNodes[1];
                    ptr->ltChild = kdNodes[0];
                    ptr->gtChild = kdNodes[2];
                } else {
                    // kdNodes[0]->tuple < kdNodes[1]->tuple; kdNodes[2]->tuple < kdNodes[1]->tuple
                    if (MergeSort<K,V>::superKeyCompare(kdNodes[0]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                        // kdNodes[0]->tuple < kdNodes[2]->tuple < kdNodes[1]->tuple
                        ptr = kdNodes[2];
                        ptr->ltChild = kdNodes[0];
                        ptr->gtChild = kdNodes[1];
                    } else {
                        // kdNodes[2]->tuple < kdNodes[0]->tuple < kdNodes[1]->tuple
                        ptr = kdNodes[0];
                        ptr->ltChild = kdNodes[2];
                        ptr->gtChild = kdNodes[1];
                    }
                }
            } else {
                // kdNodes[1]->tuple < kdNodes[0]->tuple
                if (MergeSort<K,V>::superKeyCompare(kdNodes[0]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                    // kdNodes[1]->tuple < kdNodes[0]->tuple < kdNodes[2]->tuple
                    ptr = kdNodes[0];
                    ptr->ltChild = kdNodes[1];
                    ptr->gtChild = kdNodes[2];
                } else {
                    // kdNodes[1]->tuple < kdNodes[0]->tuple; kdNodes[2]->tuple < kdNodes[0]->tuple
                    if (MergeSort<K,V>::superKeyCompare(kdNodes[1]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                        // kdNodes[1]->tuple < kdNodes[2]->tuple < kdNodes[0]->tuple
                        ptr = kdNodes[2];
                        ptr->ltChild = kdNodes[1];
                        ptr->gtChild = kdNodes[0];
                    }
                    else {
                        // kdNodes[2]->tuple < kdNodes[1]->tuple < kdNodes[0]->tuple
                        ptr = kdNodes[1];
                        ptr->ltChild = kdNodes[2];
                        ptr->gtChild = kdNodes[0];
                    }
                }
            }
            ptr->ltChild->height = ptr->gtChild->height = 1;
            ptr->height = 2;

        } else {

            // This is an illegal condition that should never occur.
            ostringstream buffer;
            buffer << "\n\n" << kdNodes.size() << " KdNode instances in rebuildSubTree1to3\n";
            throw runtime_error(buffer.str());

        }
        return ptr;
    }

    /*
     * Count the number of nodes in a subtree.     
     *
     * Calling parameter:
     * 
     * @param node (IN) the root of the subtree
     * 
     * return the number of nodes
     */
private:
    size_t countNodes(KdNode<K,V>* const node) {

        // Count this node.
        size_t count = 1;

        // Obtain counts recursively.
        if (node->ltChild != nullptr) {
            count += countNodes(node->ltChild);
        }
        if (node->gtChild != nullptr) {
            count += countNodes(node->gtChild);
        }

        return count;
    }

    /*
     * Append the nodes from a subtree into a pre-sized vector.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param kdNodes (MODIFIED) the vector
     * @param index (MODIFIED) the index where the node is stored
     */
private:
    void getSubTree(KdNode<K,V>* const node,
                    vector<KdNode<K,V>*>& kdNodes,
                    size_t& index) {

        // Append the child nodes recursively.
        if (node->ltChild != nullptr) {
            getSubTree(node->ltChild, kdNodes, index);
        }
        if (node->gtChild != nullptr) {
            getSubTree(node->gtChild, kdNodes, index);
        }

        // Initialize this node's fields, including its child
        // pointers, and append it to the kdNodes vector.
        node->height = 1;
        node->ltChild = node->gtChild = nullptr;
        kdNodes[index++] = node;
    }

#ifdef ENABLE_1TO3

    /*
     * Build a subtree that contains 1 to 3 KdNode instances
     * but skip the root node.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param nodeCount (IN) the number of nodes in the subtree
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K,V>* rebuildSubTreeSkipRoot1to3(KdNode<K,V>* const node,
                                                   size_t const nodeCount,
                                                   size_t const dim,
                                                   signed_size_t const p) {

        // Allocate a sized vector of KdNode pointers, and walk the subtree
        // to assign to each vector element a pointer to a KdNode instance.
        vector<KdNode<K,V>*> kdNodes(nodeCount);
        size_t index = 0;
        getSubTreeSkipRoot(node, kdNodes, index);
        return rebuildSubTree1to3(node, kdNodes, dim, p);
    }

    /*
     * Count the number of nodes in a subtree but skip the root node. 
     *
     * Calling parameter:
     * 
     * @param node (IN) the root of the subtree
     * @param depth (IN) the depth in the subtree (default 0)
     * 
     * return the number of nodes
     */
private:
    size_t countNodesSkipRoot(KdNode<K,V>* const node,
                              size_t const depth = 0) {

        // Initialize the count.
        size_t count = 0;

        // Obtain counts recursively.
        if (node->ltChild != nullptr) {
            count += countNodesSkipRoot(node->ltChild, depth+1);
        }
        if (node->gtChild != nullptr) {
            count += countNodesSkipRoot(node->gtChild, depth+1);
        }

        // Count this node unless it's the root.
        if (depth > 0) {
            ++count;
        }

        return count;
    }

    /*
     * Append the nodes from a subtree into a pre-sized vector
     * but skip the root node.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param kdNodes (MODIFIED) the vector
     * @param index (MODIFIED) the index where the node is stored
     * @param depth (IN) the depth in the subtree (default 0)
     */
private:
    void getSubTreeSkipRoot(KdNode<K,V>* const node,
                            vector<KdNode<K,V>*>& kdNodes,
                            size_t& index,
                            size_t const depth = 0) {

       // Append the child nodes recursively.
        if (node->ltChild != nullptr) {
            getSubTreeSkipRoot(node->ltChild, kdNodes, index, depth+1);
        }
        if (node->gtChild != nullptr) {
            getSubTreeSkipRoot(node->gtChild, kdNodes, index, depth+1);
        }

        // Initialize this node's fields, including its child
        // pointers, and append it to the kdNodes vector
        // unless it's the root node.
       if (depth > 0) {
            node->height = 1;
            node->ltChild = node->gtChild = nullptr;
            kdNodes[index++] = node;
        }
    }

#endif // ENABLE_1TO3

    /*
     * Count and append the tuples of nodes from a subtree into a pre-sized vector.
     * Note that the subtree is traversed in order; hence, the tuples are sorted.   
     *
     * Calling parameters:
     * 
     * @param coordinates (MODIFIED) a vector of key-value pairs
     * 
     * return the number of nodes
     */
public:
    size_t getSortedTree(vector<pair<vector<K>, V>>& coordinates) {

        size_t index = 0;
        return getSortedTree(KdTree<K,V>::root, coordinates, index);
    }

    /*
     * Count and append the tuples of nodes from a subtree into a pre-sized vector.
     * Note that the subtree is traversed in order; hence, the tuples are sorted.   
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param coordinates (MODIFIED) a vector of key-value pairs
     * @param index (MODIFIED) the index where the node is stored
     * 
     * return the number of nodes
     */
private:
    size_t getSortedTree(KdNode<K,V>* const node,
                         vector<pair<vector<K>, V>>& coordinates,
                         size_t& index) {

        // Initialize the count.
        size_t count = 0;

        // Obtain counts from < child nodes and copy tuples recursively.
        if (node->ltChild != nullptr) {
            count += getSortedTree(node->ltChild, coordinates, index);
        }

        // Count this node and copy its tuple and
        // the first value in its value set.
        ++count;
        for (size_t i = 0; i < coordinates[0].first.size(); ++i) {
            coordinates[index].first[i] = node->tuple[i];
        }
        coordinates[index].second = *node->values->begin();
        ++index;

        // Obtain counts from the > child nodes and copy tuples recursively.
        if (node->gtChild != nullptr) {
            count += getSortedTree(node->gtChild, coordinates, index);
        }

        return count;
    }

    /*
     * Determine whether a subtree is balanced.
     *
     * Calling parameter:
     * 
     * @param node (IN) the root of the subtree
     * 
     * return true if the subtree is balanced; otherwise, false
     */
private:
    inline static bool isBalanced(KdNode<K,V>* const node) {

        // Get and order the heights at the child nodes.
        size_t const ltHeight = getHeight(node->ltChild);
        size_t const gtHeight = getHeight(node->gtChild);
        size_t hiHeight, loHeight;
        if (ltHeight > gtHeight) {
            hiHeight = ltHeight;
            loHeight = gtHeight;
        } else {
            hiHeight = gtHeight;
            loHeight = ltHeight;
        }

#ifdef AVL_BALANCE
        // AVL balancing.
        if ( (hiHeight - loHeight) > HEIGHT_DIFF ) {
            return false;
        } else {
            return true;
        }
#else
        // Red-black balancing. Does the low child exist?
        if (loHeight == 0) {
            // No, the low child does not exist, so test the
            // balance as for a standard AVL tree wherein the
            // difference in child heights must not exceed
            // HEIGHT_DIFF, for which the default value is 1.
            if (hiHeight > HEIGHT_DIFF) {
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
#endif

    }

    /*
     * Recompute the height at the KdNode instance.
     *
     * Calling parameter:
     * 
     * @param node (IN) the KdNode instance
     * 
     * @return the recomputed height at the Kdnode instance
     */
private:
    inline static size_t computeHeight(KdNode<K,V>* const node) {

        return 1 + ( (getHeight(node->ltChild) > getHeight(node->gtChild))
                     ? getHeight(node->ltChild) : getHeight(node->gtChild) );

    }

    /*
     * Return the subtree height at the KdNode instance.
     *
     * Calling parameter:
     *
     * @param node (IN) the KdNode instance
     * 
     * @return the height at an existent node; otherwise, 0
     */
public:
    inline static size_t getHeight(KdNode<K,V>* const node) {
        if (node == nullptr) {
            return 0;
        }
        return node->height;
    }

    friend class KdTree<K,V>;
    friend class KdNode<K,V>;
};

#endif // KD_MAP_DYNAMIC_H
