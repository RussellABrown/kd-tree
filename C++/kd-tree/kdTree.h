/*
 * Copyright (c) 2015, 2021, 2023, 2024, 2025 Russell A. Brown
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

/*
 * The following compilation defines are relevant.
 *
 * -D PREALLOCATE - If defined, all instances of KdNodes and all tuple arrays are
 *                  allocated en masse within vectors instead of being allocated
 *                  individually. This decreases the time required to allocate and
 *                  deallocate the KdNode instances and tuple arrays.
 * 
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 *
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * 
 * -D MERGE_CUTOFF=n - A cutoff for switching from 1 to 2 threads to merge reference
 *                     arrays in the KdNode::mergeSort* functions (default 4096)
 * 
 * -D KNLOGN_CUTOFF=n - A cutoff for using multiple threads in buildKdTree (default 4096)
 */

#ifndef KD_TREE_H
#define KD_TREE_H

/* A cutoff for using multiple threads in buildKdTree */
#ifndef KNLOGN_CUTOFF
#define KNLOGN_CUTOFF (4096)
#endif

#include "kdTreeNode.h"

#ifdef NLOGN
#include "kdTreeNlogn.h"
#else
#if !defined(YUCAO) || defined(KD_TREE_DYNAMIC_H)
#include "kdTreeKnlogn.h"
#else
#include "kdTreeYuCao.h"
#endif
#endif

/* Forward references to all classes to support any order of compilation */
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

template <typename>
class KdTreeNlogn;

template <typename>
class KdTreeKnlogn;

template <typename>
class KdTreeYuCao;

/* The KdTree class defines the k-d tree API. */
template <typename K>
class KdTree
{
public:
  KdNode<K>* root = nullptr;
  signed_size_t numDimensions = 3;
  signed_size_t maxSubmitDepth = -1;

#ifdef KD_TREE_DYNAMIC_H

  /*
   * This constructor assigns the KdTree::root node from
   * the root node of another KdTree instance, sets the
   * root node of that other KdTree instance to nullptr
   * (which transfers that other root node to KdTree::root),
   * deletes that other KdTree instance, and sets it to
   * nullptr so that it doesn't remain a dangling pointer.
   *
   * Calling parameters:
   * 
   * numDimensions (IN) the number of dimension k of the k-d tree
   * maxSubmitDepth (IN) the maximum tree depth for creating a child thread
   * tree (MODIFIED) a pointer to a KdTree instance
   */
public:
  KdTree(signed_size_t const numDimensions,
         signed_size_t const maxSubmitDepth,
         KdTree<K>*& tree)
  {
    this->numDimensions = numDimensions;
    this->maxSubmitDepth = maxSubmitDepth;
    this->root = tree->root;
    tree->root = nullptr;
    delete tree;
    tree = nullptr;
  }

#endif

  /*
   * This is the basic constructor.
   *
   * Calling parameters:
   * 
   * numDimensions (IN) the number of dimension k of the k-d tree
   * maxSubmitDepth (IN) the maximum tree depth for creating a child thread
   */
public:
  KdTree(signed_size_t const numDimensions,
         signed_size_t const maxSubmitDepth)
  {
    this->numDimensions = numDimensions;
    this->maxSubmitDepth = maxSubmitDepth;
  }

#if defined(PREALLOCATE) && !defined(KD_TREE_DYNAMIC_H)
public:
  vector<KdNode<K>>* kdNodes;
  vector<K>* tuples;
#endif

public:
  ~KdTree() {

    // If the KdNode instances and the tuples are contained by preallocated
    // vectors, delete them; otherwise, delete the KdNode::root so that
    // the ~KdNode destructor will recursively delete all tuple arrays.
    //
    // However, do not delete the KdNode::root if KD_TREE_DYNAMIC_H
    // is defined, because the KdTreeDynamic::rebuildSubTree function
    // requires that the k-d tree persist. Instead, the ~KdTreeDynamic
    // destructor will delete KdNode::root
#ifndef KD_TREE_DYNAMIC_H
#ifdef PREALLOCATE
    delete kdNodes;
    delete tuples;
#else
    delete root;
#endif
#endif

  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   * This version of createKdTree is called from the KdTreeDynamic::rebuildSubtree function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode<K>*> whose KdNodes store the (x, y, z, w...) coordinates
   * dim - the number of dimensions (required when KD_TREE_DYNAMIC_H is defined)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
  * p - the leading dimension
   *
   * returns: a KdTree pointer
   */
public:
  static KdTree<K>* createKdTree(vector<KdNode<K>*> const& kdNodes,
                                 size_t const dim,
                                 signed_size_t const maximumSubmitDepth,
                                 signed_size_t const p)
  {

#ifdef NLOGN
    return KdTreeNlogn<K>::createKdTree(kdNodes, dim, maximumSubmitDepth, p);
#else
#if !defined(YUCAO) || defined(KD_TREE_DYNAMIC_H)
    return KdTreeKnlogn<K>::createKdTree(kdNodes, dim, maximumSubmitDepth, p);
#else
    return KdTreeYuCao<K>::createKdTree(kdNodes, dim, maximumSubmitDepth, p);
#endif
#endif

  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   *
   * Calling parameters:
   *
   * coordinates - a vector<vector<K>> that stores the (x, y, z, w...) coordinates
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * numberOfNodes - the number of nodes counted by KdNode::verifyKdTree - returned by reference
   * allocateTime, sortTime, removeTime, kdTime, verifyTime, deallocateTime, unsortTime - execution times
   * p - optional leading dimension, aka partition coordnate (default 0)
   *
   * returns: a KdTree pointer
   */
public:
  static KdTree<K>* createKdTree(vector<vector<K>> const& coordinates,
                                 signed_size_t const maximumSubmitDepth,
                                 signed_size_t& numberOfNodes,
                                 double& allocateTime,
                                 double& sortTime,
                                 double& removeTime,
                                 double& kdTime,
                                 double& verifyTime,
                                 double& deallocateTime,
                                 double& unsortTime)
  {

#ifdef NLOGN
    return KdTreeNlogn<K>::createKdTree(coordinates, maximumSubmitDepth,
                                        numberOfNodes, allocateTime,
                                        sortTime, removeTime, kdTime,
                                        verifyTime, deallocateTime,
                                        unsortTime);
#else
#if !defined(YUCAO) || defined(KD_TREE_DYNAMIC_H)
    return KdTreeKnlogn<K>::createKdTree(coordinates, maximumSubmitDepth,
                                         numberOfNodes, allocateTime,
                                         sortTime, removeTime, kdTime,
                                         verifyTime, deallocateTime,
                                         unsortTime);
#else
    return KdTreeYuCao<K>::createKdTree(coordinates, maximumSubmitDepth,
                                        numberOfNodes, allocateTime,
                                        sortTime, removeTime, kdTime,
                                        verifyTime, deallocateTime,
                                        unsortTime);
#endif
#endif

  }

  /* Determine whether the tree is empty. */
public:
  bool isEmpty() {
    return (root == nullptr);
  }
  
  /* Return the root of the tree. */
public:
  KdNode<K>* getRoot() {
    return root;
  }

  /*
   * Return the height of the tree, or 0 if the tree is empty
   * or if KD_MAP_DYNAMIC_H is not defined.
   */
public:
  size_t getHeight() {

#ifdef KD_TREE_DYNAMIC_H
    if (root == nullptr) {
      return 0;
    } else {
      return root->height;
    }
#else
    return 0;
#endif

  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list<KdNode<K>*> that is passed by reference and modified
   * queryLower - the query lower bound vector that is passed by reference and modified
   * queryUpper - the query upper bound vector that is passed by reference and modified
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * enableAll - a boolean that specifies whether to test all dimensions for insidedness
   *             and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    bool const enableAll) {

    if (root != nullptr) {
      root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth, enableAll);
    }
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list of KdNode pointers that is passed by reference and modified
   * queryLower - the query lower bound vector that is passed by reference and modified
   * queryUpper - the query upper bound vector that is passed by reference and modified
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    vector<bool> const& enable) {
    
    if (root != nullptr) {
      root->searchRegion(result, queryLower, queryUpper, maximumSubmitDepth, enable);
    }
  }

  /*
  * Verify the results of region search.
  *
  * fastRegionList - a list of KdNode pointers found by region search
  * slowRegionList - a list of KdNode pointers found by brute-force search
  */
public:
  void verifyRegionSearch(list<KdNode<K>*> fastRegionList,
                          list<KdNode<K>*> slowRegionList) {

    if (root != nullptr) {
      root->verifyRegionSearch(fastRegionList, slowRegionList);
    }
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to attempt to find
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors) {
    
    if (root != nullptr) {
      root->findNearestNeighbors(neighbors, query, numNeighbors);
    }
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                            vector<bool> const& enable) {
    
    if (root != nullptr) {
      root->findNearestNeighbors(neighbors, query, numNeighbors, enable);
    }
  }
  
  /*
   * Verify the consistency between the nearest neighbors lists found
   * by k-d tree search and by brute force.
   *
   * Calling parameters:
   *
   * neighborsFast - a list of nearest neighbors found by k-d tree search
   * neighborsSlow - a list of nearest neighbors found by brute force.
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the lists. Hence, this function is not static.
   */
public:
  void verifyNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighborsFast,
                              forward_list< pair<double, KdNode<K>*> >& neighborsSlow) const {
    
    if (root != nullptr) {
      root->verifyNearestNeighbors(neighborsFast, neighborsSlow);
    }
  }
  
  /*
   * Find M nearest neighbors to the query vector via brute force
   * and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   */
public:
  void bruteNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                             vector<K> const& query,
                             signed_size_t const numNeighbors) {
    
    if (root != nullptr) {
      root->bruteNearestNeighbors(neighbors, query, numNeighbors);
    }
  }
  
  /*
   * The verifyKdTree function walks the k-d tree and checks that the
   * children of a node are in the correct branch of that node.
   *
   * returns: a count of the number of kdNodes in the k-d tree
   */
public:
  signed_size_t verifyKdTree() {

    if (root != nullptr) {
      return root->verifyKdTree(numDimensions, maxSubmitDepth, 0, 0);
   } else {
      return 0;
    }
  }

  /*
   * Sort a list of KdNode instances by increasing distance from the query point.
   *
   * Calling parameters:
   * 
   * kdList - the list of KdNode instances
   * query - a vector that contains the query point coordinates
   * maxNodes - the maximum number of nodes to maintain on the heap
   * 
   * returns: a sorted forward list of (double, KdNode*) pairs
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  forward_list<pair<double, KdNode<K>*>> sortByDistance(list<KdNode<K>*> const& kdList,
                                                         vector<K> const& query,
                                                         signed_size_t const& maxNodes) {

    if (root != nullptr) {
      return root->sortByDistance(kdList, query, maxNodes);
    } else {
      forward_list<pair<double, KdNode<K>*>> sortedList;
      return sortedList;
    }
  }

#ifdef REVERSE_NEAREST_NEIGHBORS
  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to a nearest neighbors vector and to a reverse
   * nearest neighbors vector.
   *
   * Each element of the reverse nearest neighbors vector contains a list
   * of KdNodes to which the reference KdNode is a nearest neighbor.
   *
   * The concept of reverse nearest neighbors was first described by F. Korn and
   * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
   * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
   * Data, pages 201-212.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   */
public:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t maximumSubmitDepth) {
    
    if (root != nullptr) {
      root->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions,
                                        numNeighbors, maximumSubmitDepth);
    }
  }

  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to a nearest neighbors vector and to a reverse
   * nearest neighbors vector.
   *
   * Each element of the reverse nearest neighbors vector contains a list
   * of KdNodes to which the reference KdNode is a nearest neighbor.
   *
   * The concept of reverse nearest neighbors was first described by F. Korn and
   * S. Muthukrishnan in "Influence Sets Based on Reverse Nearest Neigbor Queries",
   * Proceedings of the 2000 ACM SIGMOD International Conference on Management of
   * Data, pages 201-212.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * enable - a vector that specifies the dimensions for which to test distance
   */
public:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> >* >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> >* >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t const maximumSubmitDepth,
                                   vector<bool> const& enable) {
    
    if (root != nullptr) {
      root->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions,
                                        numNeighbors, maximumSubmitDepth, enable);
    }
  }

  /*
   * Verify the correctness of the reverse nearest neighbors vector.
   *
   * Calling parameter:
   *
   * nn - the nearest neighbors vector
   * rnn - the reverse nearest neighbors vector
   * numberOfNodes - the number of nodes in the k-d tree
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vectors. Hence, this function is not static.
   */
public:
  void verifyReverseNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                              vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                              signed_size_t const numberOfNodes) {

    if (root != nullptr) {
      root->verifyReverseNeighbors(nn, rnn, numberOfNodes);
    }
  }

  /*
   * Calculate the mean and standard deviation of the distances and list sizes in a vector.
   *
   * Calling parameter:
   *
   * vec - a vector of lists
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vector. Hence, this function is not static.
   */
public:
  void calculateMeanStd(vector< forward_list< pair<double, KdNode<K>*> > >& vec,
                        double& meanSize,
                        double& stdSize,
                        double& meanDist,
                        double& stdDist) const {

    if (root != nullptr) {
      root->calculateMeanStd(vec, meanSize, stdSize, meanDist, stdDist);
    }
  }

  /*
   * Count the number of non-empty lists in a vector.
   *
   * Calling parameter:
   *
   * vec - a vector of lists
   *
   * Although this function does not directly access the k-d tree, it requires the persistence
   * of the k-d tree for access to the KdNodes via the vector. Hence, this function is not static.
   */
public:
  size_t nonEmptyLists(vector< forward_list< pair<double, KdNode<K>*> > >& vec) const {

    if (root != nullptr) {
      return root->nonEmptyLists(vec);
    } else {
      return 0;
    }
  }
#endif // REVERSE_NEAREST_NEIGHBORS
                                                                                               
  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameters:
   *
   * tuple - the tuple as an array
   * dim - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuple(K const* tuple,
                  signed_size_t const dim) const {
    
    if (root != nullptr) {
      root->printTuple(tuple, dim);
    }
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameter:
   *
   * tuple - the tuple as a vector
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuple(vector<K> const& tuple) const {
    
    KdNode<K>::printTuple(tuple);
  }

 /*
   * The printTuples function prints all tuples in a forward list of pairs.
   *
   * Calling parameters:
   *
   * regionList - a list of (double, KdNode*) pairs
   * maximumNumberOfNodesToPrint - the maximum number of KdNodes to print
   * numDimensions - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
public:
  void printTuples(forward_list<pair<double, KdNode<K>*>> const& regionList,
                   signed_size_t const maximumNumberOfNodesToPrint,
                   signed_size_t const numDimensions) const {
    
    if (root != nullptr) {
      root->printTuples(regionList, maximumNumberOfNodesToPrint, numDimensions);
    }
  }

  /*
   * The printKdTree function prints the k-d tree "sideways" with the root at the left.
   *
   * Calling parameters:
   *
   * dim - the number of dimensions
   */
public:
  void printKdTree(signed_size_t const dim) const {
    
    if (root != nullptr) {
      root->printKdTree(root, dim, 0);
    }
  }

}; // class KdTree

#endif // KD_TREE_H
