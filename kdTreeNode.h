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
 * The k-d tree was described by Jon Bentley in "Multidimensional Binary Search Trees
 * Used for Associative Searching", CACM 18(9): 509-517, 1975.  For k dimensions and
 * n elements of data, a balanced k-d tree is built in O(kn log n) + O((k-1)n log n)
 * time by first sorting the data in each of k dimensions, then building the k-d tree
 * in a manner that preserves the order of the k sorts while recursively partitioning
 * the data at each level of the k-d tree. No further sorting is necessary.
 */
 
#ifndef KD_TREE_NODE_H
#define KD_TREE_NODE_H

#include <algorithm>
#include <chrono>
#include <climits>
#include <exception>
#include <forward_list>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utility>
#include <vector>

using std::async;
using std::cout;
using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::endl;
using std::distance;
using std::exception;
using std::fixed;
using std::forward_list;
using std::future;
using std::launch;
using std::list;
using std::map;
using std::lock_guard;
using std::make_pair;
using std::min;
using std::mutex;
using std::numeric_limits;
using std::ostringstream;
using std::pair;
using std::ref;
using std::runtime_error;
using std::scientific;
using std::setprecision;
using std::streamsize;
using std::vector;

/* Convert microseconds to seconds for use with std::chrono */
static const double MICROSECONDS_TO_SECONDS = 1000000.;

/* This type is the signed equivalent of size_t and might be equivalent to intmax_t */
typedef streamsize signed_size_t;

#include "kdTreeMergeSort.h"
#include "kdTreeHeapSort.h"

/* A forward reference to the friend class KdTree */
template <typename>
class KdTree;

/* One node of a k-d tree where T is key type */
template <typename K>
class KdNode {
private:
  K* tuple;
  KdNode<K>* ltChild;
  KdNode<K>* gtChild;
#ifdef REVERSE_NEAREST_NEIGHBORS
  size_t index;
#endif

public:
  KdNode() {
    ltChild = gtChild = nullptr;
  }

public:
  ~KdNode() {
    delete[] tuple;
    delete ltChild;
    delete gtChild;
  }

private:
  K* getTuple() {
    return this->tuple;
  }

  /*
   * The removeDuplicates function checks the validity of the merge sort and
   * removes duplicates from a reference array.
   *
   * Calling parameters:
   *
   * reference - a K** that represents one of the reference arrays
   * i - the leading dimension for the super key
   * dim - the number of dimensions
   *
   * returns: the end index of the reference array following removal of duplicate elements
   */
private:
  inline
  static signed_size_t removeDuplicates(K** const reference,
                                        signed_size_t const i,
                                        signed_size_t const dim,
                                        signed_size_t const size) {
    
    signed_size_t end = 0;
    for (signed_size_t j = 1; j < size; ++j) {
      auto const compare = MergeSort<K>::superKeyCompare(reference[j], reference[end], i, dim);
      if (compare < 0) {
        ostringstream buffer;
        buffer << "\n\nmerge sort failure: superKeyCompare(ref[" << j << "], ref["
               << end << "], " << i << ") = " << compare << " in removeDuplicates\n";
        throw runtime_error(buffer.str());
      }
      else if (compare > 0) {
        // Keep the jth element of the reference array.
        reference[++end] = reference[j];
      } else {
        // Skip over the jth element of the reference array and delete the tuple.
        // A pointer to any tuple not skipped is copied from the reference array
        // to the kdNodes vector; hence, this tuple will be deleted by the ~KdNode
        // destructor. So there will be no need to delete tuples via their
        // pointers in this reference array or any other references array.
        delete[] reference[j];
      }
    }
    return end;
  }

  /*
   * The getKdNode function gets a KdNode from the kdNodes array and assigns the tuple field.
   * It also assigns to the index field the index of the KdNode within the kdNodes array.
   *
   * Calling parameters:
   *
   * reference - a K** that represents one of the reference arrays
   * kdNodes - a vector<KdNode*> that contains pre-allocated KdNodes
   * k - the index into both the reference and the kdNodes arrays
   *
   * returns: a KdNode pointer to the KdNode to which the tuple has been assigned
   */
private:
  inline
  static KdNode<K>* getKdNode(K** const reference,
                              vector<KdNode<K>*> const& kdNodes,
                              signed_size_t const k) {
    
    KdNode<K>* kdNode = kdNodes[k];
    kdNode->tuple = reference[k];
#ifdef REVERSE_NEAREST_NEIGHBORS
    kdNode->index = k;
#endif
    return kdNode;
  }

  /*
   * Create a permutation vector.
   *
   * Calling parameters:
   *
   * permutation - the permutation vector that is passed by reference and modified
   * numDimensions - the number of dimensions
   * numCoordinates - the number of points in the coordinates vector
   */
private:
  static
  void createPermutation(vector<signed_size_t>& permutation,
                         signed_size_t const numDimensions,
                         signed_size_t const numCoordinates) {
    
    // Determine the maximum depth of the k-d tree, which is log2(numCoordinates).
    signed_size_t size = numCoordinates;
    signed_size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // Because the partition coordinate permutes in the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // (for e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference array
    // using 0 as the leading key of the super key.
    permutation.resize(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation[i] = i % numDimensions;
    }
  }

  /*
   * The verifyKdTree function walks the k-d tree and checks that the
   * children of a node are in the correct branch of that node.
   *
   * Calling parameters:
   *
   * permutation - the permutation vector
   * dim - the number of dimensions
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   *
   * returns: a count of the number of kdNodes in the k-d tree
   */
private:
  signed_size_t verifyKdTree(vector<signed_size_t> const& permutation,
                             signed_size_t const dim,
                             signed_size_t const maximumSubmitDepth,
                             signed_size_t const depth) const {

    signed_size_t count = 1;
    if (tuple == nullptr) {
      throw runtime_error("\n\npoint is null in verifyKdTree\n");
    }

    // The partition cycles as x, y, z, w...
    signed_size_t const p = permutation[depth];

    if (ltChild != nullptr) {
      if (ltChild->tuple[p] > tuple[p]) {
        throw runtime_error("\n\nchild is > node in verifyKdTree\n");
      }
      if (MergeSort<K>::superKeyCompare(ltChild->tuple, tuple, p, dim) >= 0) {
        throw runtime_error("\n\nchild is >= node in verifyKdTree\n");
      }
    }
    if (gtChild != nullptr) {
      if (gtChild->tuple[p] < tuple[p]) {
        throw runtime_error("\n\nchild is < node in verifyKdTree\n");
      }
      if (MergeSort<K>::superKeyCompare(gtChild->tuple, tuple, p, dim) <= 0) {
        throw runtime_error("\n\nchild is <= node in verifyKdTree\n");
      }
    }

    // Verify the < branch with a child thread at as many levels of the tree as possible.
    // Create the child thread as high in the tree as possible for greater utilization.

    // Is a child thread available to verify the < branch?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so verify the < branch with the current thread.
      if (ltChild != nullptr) {
        count += ltChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
      }

      // Then verify the > branch with the current thread.
      if (gtChild != nullptr) {
        count += gtChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
      }
    }
    else {

      // Yes, so verify the < branch with a child thread. Note that a lambda
      // is required because this verifyKdTree function is not static. The use
      // of std::ref may be unnecessary in view of the [&] lambda argument specification.

      future<signed_size_t> verifyFuture;
      if (ltChild != nullptr) {
        verifyFuture =
          async(launch::async, [&] {
                                 return ltChild->verifyKdTree(ref(permutation),
                                                              dim,
                                                              maximumSubmitDepth,
                                                              depth + 1);
                               });
      }

      // And simultaneously verify the > branch with the current thread.
      signed_size_t gtCount = 0;
      if (gtChild != nullptr) {
        gtCount = gtChild->verifyKdTree(permutation, dim, maximumSubmitDepth, depth + 1);
      }

      // Wait for the child thread to finish execution.
      signed_size_t ltCount = 0;
      if (ltChild != nullptr) {
        try {
          ltCount = verifyFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for verify future in verifyKdTree\n");
        }
      }

      // Sum the counts returned by the child and current threads.
      count += ltCount + gtCount;
    }

    return count;
  }

  /*
   * The insideBounds function determines whether KdNode::tuple lies inside the
   * hyper-rectangle defined by the query lower and upper bound vectors.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *
   * return true if inside, false if outside
   */
private:
  bool insideBounds(vector<K> const& queryLower,
                    vector<K> const& queryUpper,
                    vector<bool> const& enable) const {
    
    bool inside = true;
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (enable[i] && (queryLower[i] > tuple[i] || queryUpper[i] < tuple[i])) {
        inside = false;
        break;
      }
    }
    return inside;
  }

  /*
   * The regionSearch function searches the k-d tree recursively to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list<KdNode<K>*> that is passed by reference and modified
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * permutation - vector that specifies permutation of the partition coordinate
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   * enable - a vector that specifies the dimensions on which to prune the region search
   */
private:
  void regionSearch(list<KdNode<K>*>& result,
                    vector<K> const& queryLower,
                    vector<K> const& queryUpper,
                    vector<signed_size_t> const& permutation,
                    signed_size_t const maximumSubmitDepth,
                    signed_size_t const depth,
                    vector<bool> const&  enable) {

    // The partition cycles as x, y, z, w...
    signed_size_t const p = permutation[depth];

    // If the KdNode is within the query hyper-rectangle for each of the k dimensions,
    // add the KdNode to the list of KdNodes that lie inside the hyper-cube. The
    // following loop is equivalent to the IN_REGION pseudo-Algol code proposed
    // by Jon Bentley in his CACM article.
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_front(this);
    }

    // Determine whether to search the < and > branches of the k-d tree. Although
    // the superKeyCompare function can produce a different result for the == case
    // than does comparing only the leading keys of the super-keys, that result
    // will avoid unnecessary searching of a sub-tree (at the expense of a more
    // precise super-key comparison) but the unnecessary search appears not to
    // change the outcome of this recursive regionSearch function.
    //
    // Note that if the partition dimension is not enabled, both branches are searched.
#ifdef NO_SUPER_KEY
    bool const searchLT = ltChild != nullptr && (tuple[p] >= queryLower[p] || !enable[p]);
    bool const searchGT = gtChild != nullptr && (tuple[p] <= queryUpper[p] || !enable[p]);;
#else
    bool const searchLT = ltChild != nullptr && (MergeSort<K>::superKeyCompare(tuple, queryLower.data(),
                                                                               p, queryLower.size()) >= 0
                                                || !enable[p]);
    bool const searchGT = gtChild != nullptr && (MergeSort<K>::superKeyCompare(tuple, queryUpper.data(),
                                                                               p, queryLower.size()) <= 0
                                                || !enable[p]);
#endif

    // Do both branches require searching and is a child thread available?
    if (searchLT && searchGT && maximumSubmitDepth >= 0 && depth <= maximumSubmitDepth) {

      // Yes, both branches of the tree require searching and a child thread is available,
      // so prepare to search the < branch with a child thread.
      future<void> searchFuture;

      // Search the < branch?
      if (searchLT) {
        
        // Yes, search the < branch asynchronously with a child thread.
        // A lamba is required because this regionSearch function is not
        // static. The use of std::ref may be unnecessary in view of the
        // [&] lambda argument specification.
        list<KdNode<K>*> ltResult;
        searchFuture = async(launch::async, [&] {
                                              ltChild->regionSearch(ref(ltResult),
                                                                    ref(queryLower),
                                                                    ref(queryUpper),
                                                                    ref(permutation),
                                                                    maximumSubmitDepth,
                                                                    depth + 1,
                                                                    ref(enable));
                                            });
        // Search the > branch?
        list<KdNode<K>*> gtResult;
        if (searchGT) {
          
          // Yes, search the > branch  with the master thread.
          gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Get the result of searching the < branch with the child thread.
        try {
          searchFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for search future in regionSearch\n");
        }

        // Append the results of searching the < and > branches to the result (if any) for this KdNode.
        result.splice(result.end(), ltResult);
        result.splice(result.end(), gtResult);

      } else {

        // No, don't search the < branch. Search the > branch?
        list<KdNode<K>*> gtResult;
        if (searchGT) {
          
          // Yes, search the > branch  with the master thread.
          gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Append the result of searching the > branch to the result (if any) for this KdNode.
        result.splice(result.end(), gtResult);
      }

    } else {
      
      // No, both branches do not require searching. Search the < branch with the master thread?
      if (searchLT) {
        list<KdNode<K>*> ltResult;
        ltChild->regionSearch(ltResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), ltResult);
      }

      // Search the > branch with the master thread?
      if (searchGT) {
        list<KdNode<K>*> gtResult;
        gtChild->regionSearch(gtResult, queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), gtResult);
      }

    }
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
   * size - the number of points in the coordinates vector
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
private:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    signed_size_t const size) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, queryLower.size(), size);
    
    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        auto const tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over all dimensions to obtain the resulting list of KdNodes.
    vector<bool> enable(queryLower.size(), true);
    regionSearch(result, queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
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
   * size - the number of points in the coordinates vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
private:
  void searchRegion(list<KdNode<K>*>& result,
                    vector<K>& queryLower,
                    vector<K>& queryUpper,
                    signed_size_t const maximumSubmitDepth,
                    signed_size_t const size,
                    vector<bool> const& enable) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, queryLower.size(), size);
    
    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        auto const tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over the enabled dimensions to obtain the resulting list of KdNodes.
    regionSearch(result, queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
  }

  /*
   * Walk the k-d tree recursively and append to a list each KdNode that lies inside
   * the hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list of KdNode pointers that is passed by reference and modified
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the search
   *
   * return a list of pointers to KdNodes that lie within the query hyper-rectangle.
   */
private:
  void regionBrute(list<KdNode<K>*>& result,
                   vector<K> const& queryLower,
                   vector<K> const& queryUpper,
                   vector<bool> const& enable) {

    // Append the KdNode to the list if it lies inside the query bounds.
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_front(this);
    }

    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      list<KdNode<K>*> ltResult;
      ltChild->regionBrute(ltResult, queryLower, queryUpper, enable);
      result.splice(result.end(), ltResult);
    }

    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      list<KdNode<K>*> gtResult;
      gtChild->regionBrute(gtResult, queryLower, queryUpper, enable);
      result.splice(result.end(), gtResult);
    }
  }

  /*
   * Walk the k-d tree and append to a list each KdNode that lies inside
   * the hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * result - a list of KdNode pointers that is passed by reference and modified
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   *
   * return a list of pointers to KdNodes that lie within the query hyper-rectangle.
   */
private:
  void bruteRegion(list<KdNode<K>*>& result,
                   vector<K>& queryLower,
                   vector<K>& queryUpper) {

    // Search over all dimensions.
    vector<bool> enable(queryLower.size(), true);

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        auto const tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Walk the k-d tree recursively.
    regionBrute(result, queryLower, queryUpper, enable);
  }

  /*
   * Search the k-d tree recursively for up to M nearest geometric neighbors to a query point
   * by adding them to the NearestNeighborHeap that stores up to M neighbors.  Exclude from the
   * search any branch of the tree wherein it is guaranteed that all nodes in that branch are
   * farther away than the current farthest node stored in the heap. Details of the search
   * algorithm are described by Friedman et al. in "An Algorithm for Finding Best Matches in
   * Logarithmic Expected Time", ACM Transactions on Mathematical Software, 3: 209-226, 1977.
   *
   * Calling parameters:
   *
   * heap - an instance of NearestNeighborHeap that is built relative to a query point
   * permutation - vector that specifies permutation of the partition coordinate
   * depth - depth in the k-d tree
   */
private:
  void nearestNeighbors(NearestNeighborHeap<K>& heap,
                        vector<signed_size_t> const& permutation,
                        signed_size_t const depth) {

    // The partition permutes as x, y, z, w...
    signed_size_t p = permutation[depth];

    // If query[p] < tuple[p], descend the < branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    if (heap.query[p] < tuple[p]) {
      if (ltChild != nullptr) {  // If not at the bottom of the tree, descend the < branch unconditionally.
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the > branch and then attempt to add the node to the heap.
      double const dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (gtChild != nullptr) { // If not at the bottom of the tree, descend the > branch
          gtChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current KdNode to the heap.
      }
    }
    // If query[p] > tuple[p], descend the > branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    else if (heap.query[p] > tuple[p]) {
      if (gtChild != nullptr) {  // If not at the bottom of the tree, descend the > branch unconditionally.
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the < branch and then attempt to add the node to the heap.
      double const dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (ltChild != nullptr) {
          ltChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current node to the heap.
      }
    }
    // Because query[p] == tuple[p], the probability of finding nearest neighbors is equal for both branches
    // of the tree, so descend both branches and then attempt to add the current node to the heap.
    else {
      if (ltChild != nullptr) {
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      if (gtChild != nullptr) {
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      heap.add(this);  // Attempt to add the current node to the heap.
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
   * size - the number of points in the coordinates vector
   */
private:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                            signed_size_t const size) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, query.size(), size);

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K> heap(query, numNeighbors);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and prepending it to a list.
    // Remove only the number of heap entries present.
    signed_size_t const heapDepth = heap.heapDepth();
    for (signed_size_t i = 0; i < heapDepth; ++i) {
      neighbors.push_front(heap.removeTop());
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
   * size - the number of points in the coordinates vector
   * enable - a vector that specifies the dimensions for which to test distance
   */
private:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            signed_size_t const numNeighbors,
                            signed_size_t const size,
                            vector<bool> const& enable) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, query.size(), size);

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K> heap(query, numNeighbors, enable);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and prepending it to a list.
    // Remove only the number of heap entries present.
    signed_size_t const heapDepth = heap.heapDepth();;
    for (signed_size_t i = 0; i < heapDepth; ++i) {
      neighbors.push_front(heap.removeTop());
    }
  }
  
  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified
   * query - the query vector
   * permutation - vector that specifies permutation of the partition coordinate
   * numNeighbors - the number M of nearest neighbors to attempt to find
   */
private:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            vector<signed_size_t> const& permutation,
                            signed_size_t const numNeighbors) {
    
    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K> heap(query, numNeighbors);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and prepending it to a list.
    // Remove only the number of heap entries present.
    signed_size_t const heapDepth = heap.heapDepth();;
    for (signed_size_t i = 0; i < heapDepth; ++i) {
      neighbors.push_front(heap.removeTop());
    }
  }

  /*
   * Find up to M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified
   * query - the query vector
   * permutation - vector that specifies permutation of the partition coordinate
   * numNeighbors - the number M of nearest neighbors to attempt to  find
   * enable - a vector that specifies the dimensions for which to test distance
   */
private:
  void findNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                            vector<K> const& query,
                            vector<signed_size_t> const& permutation,
                            signed_size_t const numNeighbors,
                            vector<bool> const& enable) {
    
    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K> heap(query, numNeighbors, enable);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and prepending it to a list.
    // Remove only the number of heap entries present.
    signed_size_t const heapDepth = heap.heapDepth();;
    for (signed_size_t i = 0; i < heapDepth; ++i) {
      neighbors.push_front(heap.removeTop());
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
private:
  void verifyNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighborsFast,
                              forward_list< pair<double, KdNode<K>*> >& neighborsSlow) const {
    
    auto itf1 = neighborsFast.begin();
    auto its1 = neighborsSlow.begin();
    
    // Compare the first k-d tree (fast) distance to the first brute-force (slow) distance.
    if (itf1->first != its1->first) {
      ostringstream buffer;
      buffer << "\n\nfast distance[0] = " << itf1->first << "  !=  slow distance[0] = " << its1->first << endl;
      throw runtime_error(buffer.str());
    }
    
    // Compare the first k-d tree KdNode pointer to the first brute-force KdNode pointer.
    if (itf1->second != its1->second) {
      throw runtime_error("\n\nfast KdNode*[0] != slow KdNode*[0]\n");
    }

    // Compare the remaining distances and KdNode*
    auto itf2 = itf1;
    auto its2 = its1;
    ++itf2;
    ++its2;
    signed_size_t i = 1;
    for ( ; itf2 != neighborsFast.end(); ++itf1, ++its1, ++itf2, ++its2, ++i) {
      // Ensure that the fast distances increase monotonically.
      if (itf1->first > itf2->first) {
        ostringstream buffer;
        buffer << "\n\nfast distance[" << (i-1) << "] = " << itf1->first << "  >  fast distance[" << i << "] = " << itf2->first << endl;
        throw runtime_error(buffer.str());
      }
      // Ensure that the slow distances increase monotonically.
      if (its1->first > its2->first) {
        ostringstream buffer;
        buffer << "\n\nslow distance[" << (i-1) << "] = " << its1->first << "  >  slow distance[" << i << "] = " << its2->first << endl;
        throw runtime_error(buffer.str());
      }
      // Compare the ith k-d tree distance to the ith brute-force distance.
      if (itf2->first != its2->first) {
        ostringstream buffer;
        buffer << "\n\nfast distance[" << i << "] = " << itf2->first << "  !=  slow distance[" << i << "] = " << its2->first << endl;
        throw runtime_error(buffer.str());
      }
      // Compare the ith k-d tree KdNode pointer to the ith brute-force KdNode pointer.
      if (itf2->second != its2->second) {
        ostringstream buffer;
        buffer << "\n\nfast KdNode*[" << i << "]  !=  slow KdNode*[" << i << "]\n";
        throw runtime_error(buffer.str());
      }
    }
  }
  
#ifdef REVERSE_NEAREST_NEIGHBORS
  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to nearest neighbors and reverse nearest neighbors
   * vector.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * permutation - vector that specifies permutation of the partition coordinate
   * root - the root of the k-d tree where a search for nearest neighbors must begin
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   */
private:
  void nearestNeighborsForEach(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                               vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                               vector<mutex>& mutexes,
                               vector<signed_size_t> const& permutation,
                               KdNode<K>* const root,
                               signed_size_t const numDimensions,
                               signed_size_t const numNeighbors,
                               signed_size_t const maximumSubmitDepth,
                               signed_size_t const depth) {

    // Create a query point from the KdNode's tuple, find at most the M
    // nearest neighbors to it, prepend those neighbors to a nearest
    // neighbors list, and remove the first element of that list (which is
    // the query KdNode). Use the nnList reference to improve readability
    // without copying the list.
    vector<K> const query(tuple, tuple + numDimensions);
    auto& nnList = nn[this->index];
    root->findNearestNeighbors(nnList, query, permutation, numNeighbors);
    nnList.pop_front();

    // Iterate over the remaining list of nearest neighbors and prepend
    // the query KdNode to the reverse nearest neighbors list at the
    // vector entry for the nearest neighbor. Use a std::mutex lock
    // because multiple threads may attempt to prepend simultaneously
    // to the same reverse nearest neighbors list. Because the number
    // of reverse nearest neighbors lists greatly exceeds the number
    // of threads, the probability of thread contention for the lock
    // is small, so hopefully the std::mutex lock() function is
    // efficient for the non-contention case.
    for (auto it = nnList.begin(); it != nnList.end(); ++it) {
      auto& index = it->second->index;
      {
        lock_guard<mutex> lk(mutexes[index]);
        rnn[index].push_front(make_pair(it->first, this));
      }
    }
    
    // Are child threads available to visit both branches of the tree?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so visit the < sub-tree with the master thread.
      if (ltChild != nullptr) {
        ltChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1);
      }
    
      // And then visit the > sub-tree with the master thread.
      if (gtChild != nullptr) {
        gtChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1);
      }
    } else {

      // Yes, so recursively visit the < sub-tree with a child thread.
      // A lamba is required because this nearestNeighborsForEach function
      // is not static. The use of std::ref may be unnecessary in view of
      // the [&] lambda argument specification.
      future<void> visitFuture;
      if (ltChild != nullptr) {
        visitFuture = async(launch::async, [&] {
                                             ltChild->nearestNeighborsForEach(
                                               ref(nn),
                                               ref(rnn),
                                               ref(mutexes),
                                               ref(permutation),
                                               root,
                                               numDimensions,
                                               numNeighbors,
                                               maximumSubmitDepth,
                                               depth + 1);
                                           });
      }

      // And simultaneously visit the > sub-tree with the master thread.
      if (gtChild != nullptr) {
        gtChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1);
      }

      // Wait for the child thread to finish execution.
      if (ltChild != nullptr) {
        try {
          visitFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for visit future in nearestNeighborsForEach\n");
        }
      }
    }
  }

  /*
   * Walk the k-d tree, find up to M nearest neighbors to each point in the k-d tree,
   * and add those nearest neighbors to nearest neighbors and reverse nearest neighbors
   * vector.
   *
   * Calling parameters:
   *
   * nn - the nearest neighbors vector that is passed by reference and modified
   * rnn - the reverse nearest neighbors vector that is passed by reference and modified
   * mutexes - a vector of mutexes to make individual rnn list update thread safe
   * permutation - vector that specifies permutation of the partition coordinate
   * root - the root of the k-d tree where a search for nearest neighbors must begin
   * numDimensions - the dimensionality k of the k-d tree
   * numNeighbors - the number M of nearest neighbors to attempt to find
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   * enable - a vector that specifies the dimensions for which to test distance
   */
private:
  void nearestNeighborsForEach(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                               vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                               vector<mutex>& mutexes,
                               vector<signed_size_t> const& permutation,
                               KdNode<K>* const root,
                               signed_size_t const numDimensions,
                               signed_size_t const numNeighbors,
                               signed_size_t const maximumSubmitDepth,
                               signed_size_t const depth,
                               vector<bool> const& enable) {

    // Create a query point from the KdNode's tuple, find at most the M
    // nearest neighbors to it, prepend those neighbors to a nearest
    // neighbors list, and remove the first element of that list (which is
    // the query KdNode). Use the nnList reference to improve readability
    // without copying the list.
    vector<K> const query(tuple, tuple + numDimensions);
    auto& nnList = nn[this->index];
    root->findNearestNeighbors(nnList, query, permutation, numNeighbors, enable);
    nnList.pop_front();

    // Iterate over the remaining list of nearest neighbors and prepend
    // the query KdNode to the reverse nearest neighbors list at the
    // vector entry for the nearest neighbor. Use a std::mutex lock
    // because multiple threads may attempt to prepend simultaneously
    // to the same reverse nearest neighbors list. Because the number
    // of reverse nearest neighbors lists greatly exceeds the number
    // of threads, the probability of thread contention for the lock
    // is small, so hopefully the std::mutex lock() function is
    // efficient for the non-contention case.
    for (auto it = nnList.begin(); it != nnList.end(); ++it) {
      auto& index = it->second->index;
      {
        lock_guard<mutex> lk(mutexes[index]);
        rnn[index].push_front(make_pair(it->first, this));
      }
    }
    
    // Are child threads available to visit both branches of the tree?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so visit the < sub-tree with the master thread.
      if (ltChild != nullptr) {
        ltChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1, enable);
      }
    
      // And then visit the > sub-tree with the master thread.
      if (gtChild != nullptr) {
        gtChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1, enable);
      }
    } else {

      // Yes, so recursively visit the < sub-tree with a child thread.
      // A lamba is required because this nearestNeighborsForEach function
      // is not static. The use of std::ref may be unnecessary in view of
      // the [&] lambda argument specification.
      future<void> visitFuture;
      if (ltChild != nullptr) {
        visitFuture = async(launch::async, [&] {
                                             ltChild->nearestNeighborsForEach(
                                               ref(nn),
                                               ref(rnn),
                                               ref(mutexes),
                                               ref(permutation),
                                               root,
                                               numDimensions,
                                               numNeighbors,
                                               maximumSubmitDepth,
                                               depth + 1,
                                               enable);
                                           });
      }

      // And simultaneously visit the > sub-tree with the master thread.
      if (gtChild != nullptr) {
        gtChild->nearestNeighborsForEach(nn, rnn, mutexes, permutation, root, numDimensions,
                                         numNeighbors, maximumSubmitDepth, depth + 1, enable);
      }

      // Wait for the child thread to finish execution.
      if (ltChild != nullptr) {
        try {
          visitFuture.get();
        }
        catch (exception const& e) {
          throw runtime_error("\n\ncaught exception for visit future in nearestNeighborsForEach\n");
        }
      }
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
   */
private:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t maximumSubmitDepth) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, numDimensions, nn.size());

    // Walk the k-d tree and build the nearest neighbors lists.
    nearestNeighborsForEach(nn, rnn, mutexes, permutation, this, numDimensions,
                            numNeighbors, maximumSubmitDepth, 0);
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
private:
  void findReverseNearestNeighbors(vector< forward_list< pair<double, KdNode<K>*> >* >& nn,
                                   vector< forward_list< pair<double, KdNode<K>*> >* >& rnn,
                                   vector<mutex>& mutexes,
                                   signed_size_t const numDimensions,
                                   signed_size_t const numNeighbors,
                                   signed_size_t const maximumSubmitDepth,
                                   vector<bool> const& enable) {
    
    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    vector<signed_size_t> permutation;
    createPermutation(permutation, numDimensions, nn.size());

    // Walk the k-d tree and build the nearest neighbors lists.
    nearestNeighborsForEach(nn, rnn, mutexes, permutation, this, numDimensions,
                            numNeighbors, maximumSubmitDepth, 0, enable);
  }

  /*
   * Populate the kdNodes vector.
   *
   * Calling parameter:
   * 
   * kdNodes - the kdNodes vector that is passed by reference and modified
   */
private:
  void populateKdNodes(vector<KdNode<K>*>& kdNodes) {
    kdNodes[index] = this;
    if (ltChild != nullptr) {
      ltChild->populateKdNodes(kdNodes);
    }
    if (gtChild != nullptr) {
      gtChild->populateKdNodes(kdNodes);
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
private:
  void verifyReverseNeighbors(vector< forward_list< pair<double, KdNode<K>*> > >& nn,
                              vector< forward_list< pair<double, KdNode<K>*> > >& rnn,
                              signed_size_t const numberOfNodes) {

    // Create a kdNodes vector because the createKdTree function doesn't return it.
    vector<KdNode<K>*> kdNodes(numberOfNodes);
    populateKdNodes(kdNodes);

    // Iterate through the reverse nearest neighbors vector and verify the correctness of that list.
    for (size_t i = 0; i < rnn.size(); ++i) {
      // Get the KdNode that is a nearest neighbor to all KdNodes on the list and
      // verify that it is indeed a nearest neighbor to each KdNode on the list.
      // Use the rnnList reference to improve readability without copying the list.
      auto& rnnList = rnn[i];
      for (auto rnnIt = rnnList.begin(); rnnIt != rnnList.end(); ++rnnIt) {
        // Get the nearest neighbor list for the KdNode from the nearest neighbors vector
        // and verify that the list contains the KdNode's tuple. Use the nnList reference to
        // improve readability without copying the list.
        bool match = false;
        auto& nnList = nn[rnnIt->second->index];
        for (auto nnIt = nnList.begin(); nnIt != nnList.end(); ++nnIt) {
          if (kdNodes[i] == nnIt->second) {
            match = true;
            break;
          }
        }
        if (!match) {
          throw runtime_error("\n\nnode is not a nearest neighbor\n");
        }
      }
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
private:
  void calculateMeanStd(vector< forward_list< pair<double, KdNode<K>*> > >& vec,
                        double& meanSize,
                        double& stdSize,
                        double& meanDist,
                        double& stdDist) const {

    // Count the number of vector entries that have non-empty lists
    // and sum the distances and list lengths.
    size_t count = 0;
    double sumDist = 0.0, sumDist2 = 0.0, sumSize = 0.0, sumSize2 = 0.0;
    for (size_t i = 0; i < vec.size(); ++i) {
      auto& vecList = vec[i];
      if (!vecList.empty()) {
        ++count;
        double const size = static_cast<double>(distance(vecList.begin(), vecList.end()));
        sumSize += size;
        sumSize2 += size * size;
        for (auto listIt = vecList.begin(); listIt != vecList.end(); ++listIt) {
          double const dist2 = listIt->first;
          sumDist += sqrt(dist2);
          sumDist2 += dist2;
        }
      }
    }
    double mapSize = static_cast<double>(count);
    meanSize = sumSize / mapSize;
    stdSize = sqrt((sumSize2 - (sumSize * sumSize / mapSize)) / (mapSize - 1.0));
    meanDist = sumDist / sumSize;
    stdDist = sqrt((sumDist2 - (sumDist * sumDist / sumSize)) / (sumSize - 1.0));
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
private:
  size_t nonEmptyLists(vector< forward_list< pair<double, KdNode<K>*> > >& vec) const {

    size_t count = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
      if(!vec[i].empty()) {
        ++count;
      }
    }
    return count;
  }
#endif // REVERSE_NEAREST_NEIGHBORS
                                                                                               
  /*
   * Walk the k-d tree and attempt to add each KdNode to the NearestNeighborHeap.
   *
   * Calling parameter:
   *
   * heap - an instance of NearestNeighborHeap
   */
private:
  void allNeighbors(NearestNeighborHeap<K>& heap) {

    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->allNeighbors(heap);
    }
    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->allNeighbors(heap);
    }
    // Attempt to add the current KdNode to the heap.
    heap.add(this);
  }

  /*
   * Find M nearest neighbors to the query vector via brute force and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * neighbors - the nearest neighbors list that is passed by reference and modified.
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   */
private:
  void bruteNearestNeighbors(forward_list< pair<double, KdNode<K>*> >& neighbors,
                             vector<K> const& query,
                             signed_size_t const numNeighbors) {
    
    // Create the heap, walk the k-d tree, and attempt to add each KdNode to the heap.
    NearestNeighborHeap<K> heap(query, numNeighbors);
    allNeighbors(heap);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    for (signed_size_t i = 0; i < numNeighbors; ++i) {
      neighbors.push_front(heap.removeTop());
    }
  }
  
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
private:
  void printTuple(K const* tuple,
                  signed_size_t const dim) const {
    
    cout << "(" << tuple[0] << ",";
    for (signed_size_t i = 1; i < dim - 1; ++i) cout << tuple[i] << ",";
    cout << tuple[dim - 1] << ")";
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
private:
  void printTuple(vector<K> const& tuple) const {
    
    cout << "(" << tuple[0] << ",";
    for (size_t i = 1; i < tuple.size() - 1; ++i) cout << tuple[i] << ",";
    cout << tuple[tuple.size() - 1] << ")";
  }

  /*
   * The printTuples func prints all tuples in a list.
   *
   * Calling parameters:
   *
   * regionList - a list of KdNodes returned by a region search
   * maximumNumberOfNodesToPrint - the maximum number of KdNodes to print
   * numDimensions - the number of dimensions
   *
   * Because this function does not access the k-d tree, it could be static.
   * However, calling it as a static function requires speicification of a
   * type, so calling it as a non-static function is less cumbersome.
   */
private:
  void printTuples(list< KdNode<K>* > const& regionList,
                   signed_size_t const maximumNumberOfNodesToPrint,
                   signed_size_t const numDimensions) const {
    
    if (regionList.size() != 0) {
      signed_size_t maxNodesToPrint = maximumNumberOfNodesToPrint;
      for (auto it = regionList.begin(); it != regionList.end(); ++it) {
        printTuple((*it)->getTuple(), numDimensions);
        cout << endl;
        --maxNodesToPrint;
        if (maxNodesToPrint == 0) {
          break;
        }
      }
    }
  }

  /*
   * The printKdTree function prints the k-d tree "sideways" with the root at the ltChild.
   *
   * Calling parameters:
   *
   * dim - the number of dimensions
   * depth - the depth in the k-d tree
   */
private:
  void printKdTree(signed_size_t const dim,
                   signed_size_t const depth) const {
    
    if (gtChild != nullptr) {
      gtChild->printKdTree(dim, depth + 1);
    }
    for (signed_size_t i = 0; i < depth; ++i) cout << "       ";
    printTuple(tuple, dim);
    cout << endl;
    if (ltChild != nullptr) {
      ltChild->printKdTree(dim, depth + 1);
    }
  }

  friend class KdTree<K>;
  friend class NearestNeighborHeap<K>;
}; // class KdNode

#endif // KD_TREE_NODE_H
