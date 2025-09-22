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

#ifndef KD_TREE_DYNAMIC_H
#define KD_TREE_DYNAMIC_H

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::cout;
using std::endl;
using std::ostringstream;
using std::runtime_error;
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
 * However, this kdTreeDynamic.h file must be compiled before kdTreeNode.h,
 * kdTreeKnlogn.h, and kdTreeNlogn.h because those files contain #if and #ifdef
 * directives that check whether KD_TREE_DYNAMIC_H is defined.
 */
template <typename>
class KdTreeDynamic;

template <typename>
class KdTree;

template <typename>
class KdNode;

template <typename>
class MergeSort;

template <typename>
class NearestNeighborHeap;

/*
 * The KdTreeDynamic class is derived from the KdTree class and
 * accessess KdTree::root so that an instance of KdTreeDynamic
 * may use KdTree member functions that require access to root.
 */
template <typename K>
class KdTreeDynamic : public KdTree<K>
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

#ifdef STATISTICS
        insertHistogram.resize(HISTOGRAM_SIZE);
        eraseHistogram.resize(HISTOGRAM_SIZE);
        insertBalanceHistogram.resize(HISTOGRAM_SIZE);
        eraseBalanceHistogram.resize(HISTOGRAM_SIZE);
        insertBalanceMaximum.resize(MAXIMUM_SIZE);
        eraseBalanceMaximum.resize(MAXIMUM_SIZE);
        insertBalanceMaxCnt.resize(MAXIMUM_SIZE);
        eraseBalanceMaxCnt.resize(MAXIMUM_SIZE);
#endif

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
                  KdNode<K>* const root) {

        this->maxSubmitDepth = maxSubmitDepth;
        this->cutoff = cutoff;
        KdTree<K>::root = root;

#ifdef STATISTICS
        insertHistogram.resize(HISTOGRAM_SIZE);
        eraseHistogram.resize(HISTOGRAM_SIZE);
        insertBalanceHistogram.resize(HISTOGRAM_SIZE);
        eraseBalanceHistogram.resize(HISTOGRAM_SIZE);
        insertBalanceMaximum.resize(MAXIMUM_SIZE);
        eraseBalanceMaximum.resize(MAXIMUM_SIZE);
        insertBalanceMaxCnt.resize(MAXIMUM_SIZE);
        eraseBalanceMaxCnt.resize(MAXIMUM_SIZE);
#endif

    }

    /*
     * If KD_TREE_DYNAMIC_H is defined, the ~KdTree destructor does not
     * delete KdTree::root to prevent recursive deletion of the k-d tree
     * whose persistence is required by KdTreeDynamic::rebuildSubTree, so
     * so delete KdTree::root here.
     */
public:
    ~KdTreeDynamic() {
        delete KdTree<K>::root;
    }

#ifdef STATISTICS
    // Statistics counters, etc.
public:
    size_t insertCount = 0, eraseCount = 0;
    size_t insertBalanceSum = 0, eraseBalanceSum = 0, eraseFindSum = 0;
    size_t insertSizeMaxCount = 0, eraseSizeMaxCount = 0;
    size_t insertBalanceCount = 0, eraseBalanceCount = 0, eraseFindCount = 0;
    double insertBalanceTime = 0, eraseBalanceTime = 0;
    double eraseFindTime = 0, eraseRecursiveTime = 0;
    vector<size_t> insertHistogram, eraseHistogram, copyHistogram;
    vector<size_t> insertBalanceHistogram, eraseBalanceHistogram;
    vector<size_t> insertBalanceMaximum, eraseBalanceMaximum;
    vector<size_t> insertBalanceMaxCnt, eraseBalanceMaxCnt;

    inline void clearStatistics() {
        insertBalanceSum = eraseBalanceSum = eraseFindSum = 0;
        insertBalanceTime = eraseBalanceTime = 0;
        eraseFindTime = eraseRecursiveTime = 0;
    }

    /*
     * Increment the histogram vector element indexed by the subtree size.
     *
     * Calling parameters:
     * 
     * @param size (IN) the subtree size
     * @param histogram (MODIFIED) the histogram vector
     */
private:
    inline void storeToHistogram(size_t const size,
                                 vector<size_t>& histogram) {

        // The histogram records sizes that are >= 1.
        if (size == 0) {
            return;
        }

        // If the size exceeds the size of the histogram vector,
        // increment the last element of that vector; otherwise,
        // increment the vector element indexed by size.
        if (size > histogram.size() - 1) {
            ++histogram[histogram.size() - 1];
        } else {
            ++histogram[size - 1];
        }
    }

    /*
     * Increment the histogram vector element indexed by the subtree size,
     * and record the maximum size.
     *
     * Calling parameters:
     * 
     * @param size (IN) the subtree size
     * @param histogram (MODIFIED) the histogram vector
     * @param maximumSize (MODIFIED) maximum size vector
     * @param maximumCnt (MODIFIED) maximum count vector
     */
private:
    inline void storeToHistogram(size_t const size,
                                 vector<size_t>& histogram,
                                 vector<size_t>& maximumSize,
                                 vector<size_t>& maximumCnt) {

        // The histogram records sizes that are >= 1.
        if (size == 0) {
            return;
        }

        // If the size exceeds the size of the histogram vector,
        // increment the last element of that vector; otherwise,
        // increment the vector element indexed by size.
        if (size > histogram.size() - 1) {
            ++histogram[histogram.size() - 1];
        } else {
            ++histogram[size - 1];
        }

        // Update the maximum size and count vectors.
        for (size_t i = 0; i < maximumSize.size(); ++i) {
            if (size == maximumSize[i]) {
                ++maximumCnt[i];
                break;
            }
            if (size > maximumSize[i]) {
                // Copy the maximum sizes and count downward from i.
                for (size_t j = maximumSize.size() - 1; j > i; --j) {
                    maximumSize[j] = maximumSize[j-1];
                    maximumCnt[j] = maximumCnt[j-1];
                }
                // Update maximum value at i, reset the maximum count, and quit.
                maximumSize[i] = size;
                maximumCnt[i] = 1;
                break;
            }
        }
    }
#endif

    /*
     * Search the tree for the existence of a key.
     *
     * Calling parameter:
     *
     * @param xkey (IN) the key to search for
     * 
     * @return true if the key was found; otherwise, false
     */
public:
    inline bool contains(vector<K> const& key) {

        if (KdTree<K>::root != nullptr) {
            signed_size_t dim = key.size();
            K* const tuple = const_cast<K* const>(key.data());
            return contains(KdTree<K>::root, tuple, dim, 0);
        } else {
            return false;
        }
    }
    
   /* Search the tree for the existence of a key.
    *
    * Calling parameters:
    *
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple of the node to search for
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension that is permuted cyclically
    * 
    * @return true if the key was found; otherwise, false
    */
private:
    inline bool contains(KdNode<K>* const node,
                         K* const key,
                         signed_size_t dim,
                         signed_size_t p) {

        KdNode<K>* ptr = node;
        while ( ptr != nullptr ) {
            if (MergeSort<K>::superKeyCompare(key, ptr->tuple, p, dim) < 0) {
                ptr = ptr->ltChild;
            } else if (MergeSort<K>::superKeyCompare(key, ptr->tuple, p, dim) > 0) {
                ptr = ptr->gtChild;
            } else {
                return true; // found the key
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
     * @param key (IN) the key to add to the tree
     * 
     * This function modifies KdTreeDynamic::root,
     * KdTreeDynamic::inserted, and KdTreeDynamic::changed.
     * 
     * @return true if the key was added as a new node; otherwise, false
     */
public:
    inline bool insert(vector<K> const& key) {

#ifdef STATISTICS
        insertBalanceCount = 0;
#endif

        inserted = changed = false;
        signed_size_t dim = key.size();
        if (KdTree<K>::root != nullptr) {
            K* const tuple = const_cast<K* const>(key.data());
            KdTree<K>::root = insert(KdTree<K>::root, tuple, dim, 0);

            // Has the height changed due to an insertion?
            if (changed) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(KdTree<K>::root) ) {
                    KdTree<K>::root->height = computeHeight(KdTree<K>::root);
                } else {

#ifndef STATISTICS
                KdTree<K>::root = rebuildSubTree(KdTree<K>::root, 1, dim, 0);
#else
                ++insertBalanceCount;
                auto beginTime = steady_clock::now();
                KdTree<K>::root = rebuildSubTree(KdTree<K>::root, 1, dim, 0);
                auto endTime = steady_clock::now();
                auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                insertBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                }
            }
        } else {

#ifdef STATISTICS
            ++insertCount;
#endif
            KdTree<K>::root = new KdNode<K>();
            KdTree<K>::root->height = 1;
            // The tuple member field will be deleted by the ~KdNode destructor.
            KdTree<K>::root->tuple = new K[dim];
            for (signed_size_t i = 0; i < dim; ++i) {
                KdTree<K>::root->tuple[i] = key[i];
            }
            inserted = true;
        }

#ifdef STATISTICS
        insertBalanceSum += insertBalanceCount;
        storeToHistogram(insertBalanceCount, insertHistogram);
#endif

        return inserted;
    }

    /*
     * Insert a key into the subtree.
     *
     * Calling parameters:
     *
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple of the node to insert
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension that is permuted cyclically
     * 
     * @return the root of the possibly rebalanced subtree
     */
private:
    KdNode<K>* insert(KdNode<K>* const node,
                      K* const key,
                      signed_size_t dim,
                      signed_size_t p) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // The return value because the node argument to this function is read only.
        KdNode<K>* nodePtr = node;

        // Determine which child to search recursively for an insertion point.
        if (MergeSort<K>::superKeyCompare(key, nodePtr->tuple, p, dim) < 0) {
            if (nodePtr->ltChild != nullptr) {
                nodePtr->ltChild = insert(nodePtr->ltChild, key, dim, p+1);
            } else {

#ifdef STATISTICS
                ++insertCount;
#endif
                nodePtr->ltChild = new KdNode<K>();
                nodePtr->ltChild->height = 1;
                // The tuple member field will be deleted by the ~KdNode destructor.
                nodePtr->ltChild->tuple = new K[dim];
                for (signed_size_t i = 0; i < dim; ++i) {
                    nodePtr->ltChild->tuple[i] = key[i];
                }
                inserted = changed = true;
            }
        } else if (MergeSort<K>::superKeyCompare(key, nodePtr->tuple, p, dim) > 0) {
            if (nodePtr->gtChild != nullptr) {
                nodePtr->gtChild = insert(nodePtr->gtChild, key, dim, p+1);
            } else {

#ifdef STATISTICS
            ++insertCount;
#endif
                nodePtr->gtChild = new KdNode<K>();
                nodePtr->gtChild->height = 1;
                // The tuple member field will be deleted by the ~KdNode destructor.
                nodePtr->gtChild->tuple = new K[dim];
                for (signed_size_t i = 0; i < dim; ++i) {
                    nodePtr->gtChild->tuple[i] = key[i];
                }
                inserted = changed = true;
            }
        } else {
            // The tree already contains the key, but report an
            // insertion for compatibility with kdMapDynamic.h
            // and specify that the tree height doesn't change.
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

#ifndef STATISTICS
                nodePtr = rebuildSubTree(nodePtr, 1, dim, p);
#else
                ++insertBalanceCount;
                auto beginTime = steady_clock::now();
                nodePtr = rebuildSubTree(nodePtr, 1, dim, p);
                auto endTime = steady_clock::now();
                auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                insertBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
            }
       }
       return nodePtr;
    }

    /*
     * Erase a key from the tree.
     *
     * Calling parameter:
     * 
     * @param key (IN) the key to erase from the tree
     * 
     * @return true if the key was found; otherwise, false
     */
public:
    inline bool erase(vector<K> const& key) {

#ifdef STATISTICS
        eraseBalanceCount = eraseFindCount = 0;
#endif

        erased = false;
        if (KdTree<K>::root != nullptr) {
            signed_size_t dim = key.size();
            K* const tuple = const_cast<K* const>(key.data());
            KdTree<K>::root = erase(KdTree<K>::root, tuple, dim, 0);
            
            // Has the height changed due to an erasure?
            if (KdTree<K>::root != nullptr && erased) {
                // Yes, the height has changed, so if the tree
                // remains balanced, compute the height at the
                // root node; otherwise rebuild the subtree to
                // rebalance it, which also computes the height.
                if ( isBalanced(KdTree<K>::root) ) {
                    KdTree<K>::root->height = computeHeight(KdTree<K>::root);
                } else {

#ifndef STATISTICS
                KdTree<K>::root = rebuildSubTree(KdTree<K>::root, 2, dim, 0);
#else
                ++insertBalanceCount;
                auto beginTime = steady_clock::now();
                KdTree<K>::root = rebuildSubTree(KdTree<K>::root, 2, dim, 0);
                auto endTime = steady_clock::now();
                auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                insertBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                }
            }
        }

#ifdef STATISTICS
        eraseBalanceSum += eraseBalanceCount;
        eraseFindSum += eraseFindCount;
        storeToHistogram(eraseBalanceCount, eraseHistogram);
#endif
        return erased;
    }

    /*
     * Erase a node from the subtree.
     * 
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param key (IN) the tuple of the node to erase
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension that is permuted cyclically
     * 
     * @return the root of the possibly rebalanced subtree
     */
private:
    KdNode<K>* erase(KdNode<K>* const node,
                     K* const key,
                     signed_size_t dim,
                     signed_size_t p) {

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

        // The return value because the node argument to this function is read only.
        KdNode<K>* nodePtr = node;

        // Determine which child to search recursively for a deletion point.
        if (MergeSort<K>::superKeyCompare(key, nodePtr->tuple, p, dim) < 0) {
            if (nodePtr->ltChild != nullptr) {
                nodePtr->ltChild = erase(nodePtr->ltChild, key, dim, p+1);

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
#ifndef STATISTICS
                        nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
#else
                        ++eraseBalanceCount;
                        auto beginTime = steady_clock::now();
                        nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
                        auto endTime = steady_clock::now();
                        auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                    }
                }
            } else {
                // The tree does not contain the node.
                erased = false;
            }
        } else if (MergeSort<K>::superKeyCompare(key, nodePtr->tuple, p, dim) > 0) {
            if (nodePtr->gtChild != nullptr) {
                nodePtr->gtChild = erase(nodePtr->gtChild, key, dim, p+1);
             
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
#ifndef STATISTICS
                        nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
#else
                        ++eraseBalanceCount;
                        auto beginTime = steady_clock::now();
                        nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
                        auto endTime = steady_clock::now();
                        auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                    }
                }
            } else {
                // The tree does not contain the node.
                erased = false;
            }
        } else {
            // Found the key, so a node will be erased.
            erased = true;

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
                    KdNode<K>* const tempPtr = nodePtr;
                    nodePtr = rebuildSubTree1to3(nodePtr->ltChild, nodeCount, 2, dim, p);
                    if (tempPtr == nullptr) {
                        ostringstream buffer;
                        buffer << "< child null in erase 1 to 3" << endl;
                        throw runtime_error(buffer.str());
                    } else {
#ifdef STATISTICS
                        ++eraseCount;
                        ++eraseBalanceCount;
#endif
                        tempPtr->ltChild = nullptr; // Prevent recursive deletion.
                        delete tempPtr;
                    }
                } else

#endif // ENABLE_1TO3

                {
                    // No, the < child subtree contains > 3 nodes. So, find
                    // the immediate predecessor node, copy the tuple from that
                    // predecessor node to the one-child node, delete the
                    // predecessor node, recompute the height along the path
                    // back to the < child (including that child), and then
                    // recompute the height at the one-child node.
                    KdNode<K>* predecessor = nodePtr->ltChild;
#ifndef STATISTICS
                    predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                    for (signed_size_t i = 0; i < dim; ++i) {
                        nodePtr->tuple[i] = predecessor->tuple[i];
                    }
                    nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, dim, p+1);
#else
                    ++eraseFindCount;
                    auto beginTime = steady_clock::now();
                    predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                    auto endTime = steady_clock::now();
                    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                    eraseFindTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
                    for (signed_size_t i = 0; i < dim; ++i) {
                        nodePtr->tuple[i] = predecessor->tuple[i];
                    }
                    beginTime = steady_clock::now();
                    nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, dim, p+1);
                    endTime = steady_clock::now();
                    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                    eraseRecursiveTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                    // Compute the height at his node because erase computes
                    // the height at the < child but not at this node. Note that,
                    // because the > child is null, deletion of a node from the
                    // < child subtree can only improve the balance. So there is no
                    // need to check the balance. However, check the balance and
                    // report an error if the subtree rooted at this node is unbalanced.
                    nodePtr->height = computeHeight(nodePtr);
                    if ( !isBalanced(nodePtr) ) {
                        ostringstream buffer;
                        buffer << "subtree is unbalanced after deletion of node from < child" << endl;
                        throw runtime_error(buffer.str());
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
                    KdNode<K>* const tempPtr = nodePtr;
                    nodePtr = rebuildSubTree1to3(nodePtr->gtChild, nodeCount, 2, dim, p);
                    if (tempPtr == nullptr) {
                        ostringstream buffer;
                        buffer << "> child null in erase 1to3" << endl;
                        throw runtime_error(buffer.str());
                    } else {
#ifdef STATISTICS
                        ++eraseCount;
                        ++eraseBalanceCount;
#endif
                        tempPtr->gtChild = nullptr; // Prevent recursive deletion.
                        delete tempPtr;
                    }
                } else

#endif // ENABLE_1TO3

                {
                    // No, the > child subtree contains > 3 nodes. So, find
                    // the immediate successor node, copy the tuple from that
                    // successor node to the one-child node, delete the
                    // successor node, recompute the height along the path
                    // back to the > child (including that child), and then
                    // recompute the height at the one-child node.
                    KdNode<K>* successor = nodePtr->gtChild;
#ifndef STATISTICS
                    successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                    for (signed_size_t i = 0; i < dim; ++i) {
                        nodePtr->tuple[i] = successor->tuple[i];
                    }
                    nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, dim, p+1);
#else
                    ++eraseFindCount;
                    auto beginTime = steady_clock::now();
                    successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                    auto endTime = steady_clock::now();
                    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                    eraseFindTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
                    for (signed_size_t i = 0; i < dim; ++i) {
                        nodePtr->tuple[i] = successor->tuple[i];
                    }
                    beginTime = steady_clock::now();
                    nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, dim, p+1);
                    endTime = steady_clock::now();
                    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                    eraseRecursiveTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                    // Compute the height at his node because erase computes
                    // the height at the > child but not at this node. Note that,
                    // because the < child is null, deletion of a node from the
                    // > child subtree can only improve the balance. So there is no
                    // need to check the balance. However, check the balance and
                    // report an error if the subtree rooted at this node is unbalanced.
                    nodePtr->height = computeHeight(nodePtr);
                    if ( !isBalanced(nodePtr) ) {
                        ostringstream buffer;
                        buffer << "subtree is unbalanced after deletion of node from > child" << endl;
                        throw runtime_error(buffer.str());
                    }
                }
            }
            // If the node has no children, delete the node.
            else if (nodePtr->ltChild == nullptr && nodePtr->gtChild == nullptr) {
                KdNode<K>* const tempPtr = nodePtr;
                nodePtr = nullptr;
                if (tempPtr == nullptr) {
                    ostringstream buffer;
                    buffer << "= = child null in erase" << endl;
                    throw runtime_error(buffer.str());
                } else {

#ifdef STATISTICS
                ++eraseCount;
#endif
                    delete tempPtr;
                }
            } else { 
                // The node has two children, so does the subtree rooted
                // at this node have <= 3 nodes, excluding this root node?
                //
                // Checking the height prior to checking the countNodesSkipRoot result
                // avoids a time-consuming call of countNodesSkipRoot for a large subtree.

#ifdef ENABLE_1TO3

                size_t nodeCount;
                if (nodePtr->height <= 3
                    && (nodeCount = countNodesSkipRoot(nodePtr)) <= 3)
                {
                    // Yes, the subtree rooted at this node contains <= 3 nodes,
                    // excluding this root node, so rebuild the subtree excluding
                    // this root node and delete this root node. This approach
                    // avoids the need to find a predecessor or successor node.
                    //
                    // Note that rebuildSubTreeSkipRoot1to3 computes the height
                    // of the rebuilt subtree.
                    KdNode<K>* const tempPtr = nodePtr;
                    nodePtr = rebuildSubTreeSkipRoot1to3(nodePtr, nodeCount, 2, dim, p);
                    if (tempPtr == nullptr) {
                        ostringstream buffer;
                        buffer << "!= != child null in erase 1to3" << endl;
                        throw runtime_error(buffer.str());
                    } else {
#ifdef STATISTICS
                        ++eraseCount;
                        ++eraseBalanceCount;
#endif
                        tempPtr->ltChild = tempPtr->gtChild = nullptr; // Prevent recursive deletion.
                        delete tempPtr;
                    }
                } else

#endif // ENABLE_1TO3

                {
                    // No, the subtree rooted at this node contains > 3 nodes,
                    // excluding this root node, so replace this node by either
                    // its predecessor or its successor.
                    //
                    // If ENABLE_PREFERRED_TEST is defined, select the replacement
                    // node from the taller of the child subtrees.

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                    cout << "< height = " << getHeight(nodePtr->ltChild)
                        << "  > height = " << getHeight(nodePtr->gtChild)
                        << "  p = " << p << endl << endl;
#endif

#ifdef ENABLE_PREFERRED_TEST

                    if ( getHeight(nodePtr->ltChild) >= getHeight(nodePtr->gtChild) )
                    {

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                        cout << "erase: must find predecessor for ";
                        KdTree<K>::printTuple(nodePtr->tuple, dim);
                        cout << endl << endl;
#endif

                        // Find the node with the largest super-key in the
                        // subtree rooted at the < child, which is the
                        // predecessor node. Copy the predecessor node's tuple
                        // to the two-child node, delete the predecessor node,
                        // recompute the heights along the path from the
                        // predecessor node to (but excluding) the two-child node,
                        // and then recompute the height at the two-child node.
                        KdNode<K>* predecessor = nodePtr->ltChild;
#ifndef STATISTICS
                        predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                        for (signed_size_t i = 0; i < dim; ++i) {
                            nodePtr->tuple[i] = predecessor->tuple[i];
                        }
                        nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, dim, p+1);
#else
                        ++eraseFindCount;
                        auto beginTime = steady_clock::now();
                        predecessor = findPredecessor(nodePtr->ltChild, predecessor, dim, p, p+1);
                        auto endTime = steady_clock::now();
                        auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseFindTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
                        for (signed_size_t i = 0; i < dim; ++i) {
                            nodePtr->tuple[i] = predecessor->tuple[i];
                        }
                        beginTime = steady_clock::now();
                        nodePtr->ltChild = erase(nodePtr->ltChild, nodePtr->tuple, dim, p+1);
                        endTime = steady_clock::now();
                        duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseRecursiveTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                        nodePtr->height = computeHeight(nodePtr);
                        if ( !isBalanced(nodePtr) ) {
#ifndef STATISTICS
                            nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
#else
                            ++eraseBalanceCount;
                            beginTime = steady_clock::now();
                            nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
                            endTime = steady_clock::now();
                            duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                            eraseBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                        }

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                        cout << "erase: predecessor = ";
                        KdTree<K>::printTuple(nodePtr->tuple, dim);
                        cout << endl << endl;
#endif
                    } else

#endif // ENABLE_PREFERRED_TEST

                    {

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                        cout << "erase: must find successor for ";
                        KdTree<K>::printTuple(nodePtr->tuple, dim);
                        cout << endl << endl;
#endif

                        // Find the node with the smallest super-key in the
                        // subtree rooted at the > child, which is the
                        // successor node. Copy the successor node's tuple
                        // to the two-child node, delete the successor node,
                        // recompute the heights along the path from the
                        // successor node to (but excluding) the two-child node,
                        // and then recompute the height at the two-child node.
                        KdNode<K>* successor = nodePtr->gtChild;
#ifndef STATISTICS
                        successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                        for (signed_size_t i = 0; i < dim; ++i) {
                            nodePtr->tuple[i] = successor->tuple[i];
                        }
                        nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, dim, p+1);
#else
                        ++eraseFindCount;
                        auto beginTime = steady_clock::now();
                        successor = findSuccessor(nodePtr->gtChild, successor, dim, p, p+1);
                        auto endTime = steady_clock::now();
                        auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseFindTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
                        for (signed_size_t i = 0; i < dim; ++i) {
                            nodePtr->tuple[i] = successor->tuple[i];
                        }
                        beginTime = steady_clock::now();
                        nodePtr->gtChild = erase(nodePtr->gtChild, nodePtr->tuple, dim, p+1);
                        endTime = steady_clock::now();
                        duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                        eraseRecursiveTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                        nodePtr->height = computeHeight(nodePtr);
                        if ( !isBalanced(nodePtr) ) {
#ifndef STATISTICS
                            nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
#else
                            ++eraseBalanceCount;
                            beginTime = steady_clock::now();
                            nodePtr = rebuildSubTree(nodePtr, 2, dim, p);
                            endTime = steady_clock::now();
                            duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
                            eraseBalanceTime += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
#endif
                        }

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                        cout << "erase: successor = ";
                        KdTree<K>::printTuple(nodePtr->tuple, dim);
                        cout << endl << endl;
#endif
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
    KdNode<K>* findPredecessor(KdNode<K>* const node,
                               KdNode<K>* const predecessor,
                               signed_size_t const dim,
                               signed_size_t const p0,
                               signed_size_t p) {

        // The return value because the predecessor argument to this function is read only.
        KdNode<K>* pred = predecessor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "findPredecessor: dim = " << dim << "  p = " << p << "  node = ";
        KdTree<K>::printTuple(node->tuple, dim);
        cout << endl << endl;
#endif

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
            // then follow both children.
            pred = checkPredecessor(node, pred, dim, p0);
            if (node->ltChild != nullptr) {
                pred = findPredecessor(node->ltChild, pred, dim, p0, p+1);

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                cout << "findPredecessor: predecessor = ";
                KdTree<K>::printTuple(pred->tuple, dim);
                cout << endl << endl;
#endif
            }
            if (node->gtChild != nullptr) {
                pred = findPredecessor(node->gtChild, pred, dim, p0, p+1);

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
            cout << "findPredecessor: predecessor = ";
            KdTree<K>::printTuple(pred->tuple, dim);
            cout << endl << endl;
#endif
            }
        }
        return pred;
    }

    /*
     * Check a potential predecessor node to see whether it is
     * a better predecessor than the current predecessor node.
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
    KdNode<K>* checkPredecessor(KdNode<K>* const node,
                                KdNode<K>* const predecessor,
                                signed_size_t const dim,
                                signed_size_t const p0) {

        // The return value because the predecessor argument to this function is read only.
        KdNode<K>* pred = predecessor;

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "checkPredecessor: p0 = " << p0
             << "  dim = " << dim << "  node = ";
        KdTree<K>::printTuple(node->tuple, dim);
        cout << "  predecessor = ";
        KdTree<K>::printTuple(pred->tuple, dim);
        cout << endl << endl;
#endif

        // Update the potential precedessor node if necessary.
        if (MergeSort<K>::superKeyCompare(node->tuple, pred->tuple, p0, dim) > 0) {
            pred = node;
        }

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "checkPredecessor: predecessor = ";
        KdTree<K>::printTuple(pred->tuple, dim);
        cout << endl << endl;
#endif
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
    KdNode<K>* findSuccessor(KdNode<K>* const node,
                             KdNode<K>* const successor,
                             signed_size_t const dim,
                             signed_size_t const p0,
                             signed_size_t p) {

        // The return value because the successor argument to this function is read only.
        KdNode<K>* succ = successor;

        // Permute the most significant dimension p cyclically using
        // a fast alternative to the modulus operator for p <= dim.
        p = (p < dim) ? p : 0;

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "findSuccessor: dim = " << dim << "  p = " << p << "  node = ";
        KdTree<K>::printTuple(node->tuple, dim);
        cout << endl << endl;
#endif

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

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                cout << "findSuccessor: successor = ";
                KdTree<K>::printTuple(succ->tuple, dim);
                cout << endl << endl;
#endif
            }
            if (node->gtChild != nullptr) {
                succ = findSuccessor(node->gtChild, succ, dim, p0, p+1);

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
                cout << "findSuccessor: successor = ";
                KdTree<K>::printTuple(succ->tuple, dim);
                cout << endl << endl;
#endif
            }
        }
        return succ;
    }

    /*
     * Check a potential successor node to see whether it is
     * a better successor than the current successor node.
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
    KdNode<K>* checkSuccessor(KdNode<K>* const node,
                              KdNode<K>* const successor,
                              signed_size_t const dim,
                              signed_size_t p0) {

        // The return value because the successor argument to this function is read only.
        KdNode<K>* succ = successor;

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "checkSuccessor: p0 = " << p0
             << "  dim = " << dim << "  node = ";
        KdTree<K>::printTuple(node->tuple, dim);
        cout << "  successor = ";
        KdTree<K>::printTuple(succ->tuple, dim);
        cout << endl << endl;
#endif

        // Update the potential successor node if necessary.
        if (MergeSort<K>::superKeyCompare(node->tuple, succ->tuple, p0, dim) < 0) {
            succ = node;
        }

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "checkSuccessor: successor = ";
        KdTree<K>::printTuple(succ->tuple, dim);
        cout << endl << endl;
#endif
        return succ;
    }

    /*
     * Rebuild a subtree to rebalance it.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param histogram (IN) specify a histogram to record statistics
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the balanced subtree
     */
private:
    inline KdNode<K>* rebuildSubTree(KdNode<K>* const node,
                                     size_t const histogram,
                                     size_t const dim,
                                     signed_size_t const p) {

        // Walk the subtree to obtain a vector of pointers to KdNodes.
        // Count the nodes in the subtree, allocate a sized vector of KdNode
        // pointers, and assign to each vector element a pointer to each KdNode.
        size_t count = countNodes(node);
        vector<KdNode<K>*> kdNodes(count);
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
             return rebuildSubTree1to3(node, kdNodes, histogram, dim, p);
        } else

#endif
        
        {
            // No, the subtree contains more than 3 KdNodes,
            // so rebuild the subtree via KdTree::createKdTree.

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
            cout << "tree prior to rebalancing subtree (p = " << p << "):" << endl << endl;
            KdTree<K>::printKdTree(dim);
            cout << endl << endl;

            cout << "must rebalance subtree (p = " << p << "):" << endl << endl;
            node->printKdTree(node, dim, 0);
            cout << endl << endl;
#endif

            // Check the count and kdNodes.size() for consistency.
            if (count != kdNodes.size()) {
                ostringstream buffer;
                buffer << "count = " << count << "  size = " << kdNodes.size() << endl;
                throw runtime_error(buffer.str());
            }

#ifdef STATISTICS
            if (histogram == 1) {
                storeToHistogram(count,
                                 insertBalanceHistogram,
                                 insertBalanceMaximum,
                                 insertBalanceMaxCnt);
            } else if (histogram == 2) {
                storeToHistogram(count,
                                 eraseBalanceHistogram,
                                 eraseBalanceMaximum,
                                 eraseBalanceMaxCnt);
            } else {
                ostringstream buffer;
                buffer << "\n\nunsupported histogram = " << histogram << " in rebuildSubTree\n";
                throw runtime_error(buffer.str());
            }
#endif

            // Call KdTree::createKdTree to rebuild the subtree, which
            // invalidates the node argument to this rebuildSubTree function.
            //
            // Do not use mutiple threads to build the subtree, even
            // though multiple threads area available, if the size of the
            // subtree is too small to justify spawning child threads.
            signed_size_t numNodes;
            double allocateTime, sortTime, removeTime, kdTime,
                verifyTime, deallocateTime, unsortTime;

            KdTree<K>* subTree;
            if (count < cutoff) {
                subTree =
                    KdTree<K>::createKdTree(kdNodes, dim, -1, numNodes,
                                            allocateTime, sortTime, removeTime, kdTime,
                                            verifyTime, deallocateTime, unsortTime, p);
            } else {
                subTree =
                    KdTree<K>::createKdTree(kdNodes, dim, maxSubmitDepth, numNodes,
                                            allocateTime, sortTime, removeTime, kdTime,
                                            verifyTime, deallocateTime, unsortTime, p);
                }

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
            cout << "rebalanced subtree (p = " << p << "):" << endl << endl;
            subTree->printKdTree(dim);
            cout << endl << endl;
#endif

            // Delete the KdTree instance but not its root
            // (see the ~KdTree destructor) so that the 
            // rebuilt subtree will be added to the k-d tree.
            KdNode<K>* const subTreeRoot = subTree->root;
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
     * @param histogram (IN) specify a histogram to record statistics
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K>* rebuildSubTree1to3(KdNode<K>* const node,
                                         size_t const nodeCount,
                                         size_t const histogram,
                                         size_t const dim,
                                         signed_size_t const p) {

        // Allocate a sized vector of KdNode pointers, and walk the subtree
        // to assign to each vector element a pointer to a KdNode instance.
        vector<KdNode<K>*> kdNodes(nodeCount);
        size_t index = 0;
        getSubTree(node, kdNodes, index);
        return rebuildSubTree1to3(node, kdNodes, histogram, dim, p);
    }

    /*
     * Build a subtree that contains 1 to 3 KdNode instances.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param kdNodes (IN) a vector of kdNode pointers
     * @param histogram (IN) specify a histogram to record statistics
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K>* rebuildSubTree1to3(KdNode<K>* const node,
                                         vector<KdNode<K>*> const& kdNodes,
                                         size_t const histogram,
                                         size_t const dim,
                                         signed_size_t const p) {

        // The return value because the node argument to this function is read only.
        KdNode<K>* ptr = node;

        // Prevent the compiler from issuing a warning when STATISTICS is undefined.
        size_t const foo = histogram;

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "rebuildSubTree1to3: node = ";
        KdTree<K>::printTuple(node->tuple, dim);
        cout << "  subtree size = " << kdNodes.size() << endl << endl;
#endif

        if (kdNodes.size() == 1) {

            // The subtree contains 1 node, so it is the root of the rebuilt subtree.
            ptr = kdNodes[0];
            ptr->height = 1;

        } else if (kdNodes.size() == 2) {

            // The subtree contains 2 nodes, so the first node is the root of the rebuilt
            // subtree. Compare super-keys to determine whether second node is the < child
            // or the > child.
            ptr = kdNodes[0];
            ptr->height = 2;
            if (MergeSort<K>::superKeyCompare(kdNodes[0]->tuple, kdNodes[1]->tuple, p, dim) > 0) {
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
            if (MergeSort<K>::superKeyCompare(kdNodes[0]->tuple, kdNodes[1]->tuple, p, dim) < 0) {
                // kdNodes[0]->tuple < kdNodes[1]->tuple
                if (MergeSort<K>::superKeyCompare(kdNodes[1]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                    // kdNodes[0]->tuple < kdNodes[1]->tuple < kdNodes[2]->tuple
                    ptr = kdNodes[1];
                    ptr->ltChild = kdNodes[0];
                    ptr->gtChild = kdNodes[2];
                } else {
                    // kdNodes[0]->tuple < kdNodes[1]->tuple; kdNodes[2]->tuple < kdNodes[1]->tuple
                    if (MergeSort<K>::superKeyCompare(kdNodes[0]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
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
                if (MergeSort<K>::superKeyCompare(kdNodes[0]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
                    // kdNodes[1]->tuple < kdNodes[0]->tuple < kdNodes[2]->tuple
                    ptr = kdNodes[0];
                    ptr->ltChild = kdNodes[1];
                    ptr->gtChild = kdNodes[2];
                } else {
                    // kdNodes[1]->tuple < kdNodes[0]->tuple; kdNodes[2]->tuple < kdNodes[0]->tuple
                    if (MergeSort<K>::superKeyCompare(kdNodes[1]->tuple, kdNodes[2]->tuple, p, dim) < 0) {
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

#ifdef STATISTICS
        if (histogram == 1) {
            storeToHistogram(kdNodes.size(),
                             insertBalanceHistogram,
                             insertBalanceMaximum,
                             insertBalanceMaxCnt);
        } else if (histogram == 2) {
            storeToHistogram(kdNodes.size(),
                             eraseBalanceHistogram,
                             eraseBalanceMaximum,
                             eraseBalanceMaxCnt);
    } else {
            ostringstream buffer;
            buffer << "\n\nunsupported histogram = " << histogram << " in rebuildSubTree1to3\n";
            throw runtime_error(buffer.str());
        }
#endif

#if defined(DEBUG_PRINT) && defined(EXTRA_PRINT)
        cout << "rebuildSubTree1to3: return node = ";
        KdTree<K>::printTuple(ptr->tuple, dim);
        cout << endl << endl;
#endif

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
    size_t countNodes(KdNode<K>* const node) {

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
    void getSubTree(KdNode<K>* const node,
                    vector<KdNode<K>*>& kdNodes,
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

    /*
     * Build a subtree that contains 1 to 3 KdNode instances
     * but skip the root node.     
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param nodeCount (IN) the number of nodes in the subtree
     * @param histogram (IN) specify a histogram to record statistics
     * @param dim (IN) the number of dimensions
     * @param p (IN) the leading dimension
     * 
     * return the root node of the subtree
     */
private:
    inline KdNode<K>* rebuildSubTreeSkipRoot1to3(KdNode<K>* const node,
                                                 size_t const nodeCount,
                                                 size_t const histogram,
                                                 size_t const dim,
                                                 signed_size_t const p) {

        // Allocate a sized vector of KdNode pointers, and walk the subtree
        // to assign to each vector element a pointer to a KdNode instance.
        vector<KdNode<K>*> kdNodes(nodeCount);
        size_t index = 0;
        getSubTreeSkipRoot(node, kdNodes, index);
        return rebuildSubTree1to3(node, kdNodes, histogram, dim, p);
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
    size_t countNodesSkipRoot(KdNode<K>* const node,
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
    void getSubTreeSkipRoot(KdNode<K>* const node,
                            vector<KdNode<K>*>& kdNodes,
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

    /*
     * Count and append the tuples of nodes from a subtree into a pre-sized vector.
     * Note that the subtree is traversed in order; hence, the tuples are sorted.   
     *
     * Calling parameters:
     * 
     * @param coordinates (MODIFIED) the vector
     * 
     * return the number of nodes
     */
public:
    size_t getSortedTree(vector<vector<K>>& coordinates) {

        size_t index = 0;
        return getSortedTree(KdTree<K>::root, coordinates, index);
    }

    /*
     * Count and append the tuples of nodes from a subtree into a pre-sized vector.
     * Note that the subtree is traversed in order; hence, the tuples are sorted.   
     *
     * Calling parameters:
     * 
     * @param node (IN) the root of the subtree
     * @param coordinates (MODIFIED) the vector
     * @param index (MODIFIED) the index where the node is stored
     * 
     * return the number of nodes
     */
private:
    size_t getSortedTree(KdNode<K>* const node,
                         vector<vector<K>>& coordinates,
                         size_t& index) {

        // Initialize the count.
        size_t count = 0;

        // Obtain counts from < child nodes and copy tuples recursively.
        if (node->ltChild != nullptr) {
            count += getSortedTree(node->ltChild, coordinates, index);
        }

        // Count this node and copy its tuple.
        ++count;
        for (size_t i = 0; i < coordinates[0].size(); ++i) {
            coordinates[index][i] = node->tuple[i];
        }
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
    inline static bool isBalanced(KdNode<K>* const node) {

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
    inline static size_t computeHeight(KdNode<K>* const node) {

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
    inline static size_t getHeight(KdNode<K>* const node) {
        if (node == nullptr) {
            return 0;
        }
        return node->height;
    }

    friend class KdTree<K>;
    friend class KdNode<K>;
};

#endif // KD_TREE_DYNAMIC_H
