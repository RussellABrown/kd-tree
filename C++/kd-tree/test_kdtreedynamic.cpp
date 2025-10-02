/*
 * Copyright (c) 2025 Russell A. Brown
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
 * Test program for kdTreeDynamic.h, kdTreeKnlogn.h and kdTreeNlogn.h
 *
 * Compile via: g++ -O3 -std=c++20 -pthread -W test_kdtreedynamic.cpp
 *
 * Optional compilation defines are as follows.
 * 
 * -D NLOGN - Select the O(n log n) algorithm instead of the O(kn log n) algorithm.
 * 
 * -D ENABLE_PREFERRED_TEST - Enable comparing the heights of a deleted 2-child node's
 *                            child subtrees to select a preferred replacement node.
 * 
 * -D ENABLE_1TO3 - Enable curtailing recursive deletion when a subtree contains <= 3 nodes.
 * 
 * -D AVL_BALANCE - If defined, KdTreeDynamic::isBalanced checks for AVL balancing;
 *                  otherwise, KdTreeDynamic::isBalanced checks for red-black balancing.
 * 
 * -D HEIGHT_DIFF=n - For red-black balancing, the maximum allowed height difference
 *                    between the < and > subtrees of a node when one subtree is empty;
 *                    for AVL balancing, the maximum allowed height difference between
 *                    the < and > subtrees of a node (default 1)
 * 
 * -D MULTI_THREAD_CUTOFF = A cutoff for multi-threaded execution of KdTree::createKdTree (default 16384)
 * 
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 *
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the KdNode::mergeSort* functions (default 15)
 * 
 * -D MERGE_CUTOFF=n - A cutoff for using multiple threads in MergeSort::mergeSort* (default 4096)
 * 
 * -D REVERSE_NEAREST_NEIGHBORS - Enable the construction of a reverse nearest neighbors
 *                                list in response to the -r command-line option.
 * 
 * The following compilation defines apply only to the O(n log n) algorithm.
 *
 * -D MEDIAN_OF_MEDIANS_CUTOFF=n - A cutoff for switching from median of medians to insertion sort
 *                                 in KdNode::partition (default 15)
 * 
 * -D MEDIAN_CUTOFF=n - A cutoff for switching from to 2 threads to calculate the median
 *                      in KdNode::partition (default 16384)
 * 
 * -D INDEX_CUTOFF=n - A cutoff for switching from to 2 threads to find the index of
 *                     the calculated median in KdNode::partition (default 16384)
 * 
 * -D BIDIRECTIONAL_PARTITION - Partition an array about the median of medians proceeding
 *                              from both ends of the array instead of only the beginning.
 * 
 * -D NLOGN_CUTOFF=n - A cutoff for using multiple threads in buildKdTree (default 4096)
 * 
 * The following compilation define applies only to the O(kn log n) algorithm.
 * 
 * -D KNLOGN_CUTOFF=n - A cutoff for using multiple threads in buildKdTree (default 4096)
 * 
 * 
 * Usage:
 *
 * test_kdtree [-i I] [-n N] [-m M] [-d D] [-t T] [-s S] [-p P] [-b] [-g] [-v] [-f] [-w] [-r]
 *
 * where the command-line options are interpreted as follows.
 * 
 * -i The number I of iterations of k-d tree insertion, search, and deletion (default 1)
 *
 * -n The number N of randomly generated points used to build the k-d tree (default 262144)
 *
 * -m The maximum number M of nearest neighbors added to a priority queue
 *    when searching the k-d tree for nearest neighbors (default 5)
 *
 * -d The number of dimensions D (aka k) of the k-d tree (default 3)
 *
 * -t The number of threads T used to build and search the k-d tree (default 1)
 *
 * -s The search divisor S used for region search (default 10)
 *
 * -p The maximum number P of nodes to report when reporting nearest neighbor results (default 5)
 * 
 * -b Build a balanced k-d tree for comparison to the dynamic k-d tree (default off)
 * 
 * -c The multi-thread cutoff below which only single-threaded execution occurs (default 65536)
 *
 * -g Find nearest neighbors to each point
 *
 * -v Verify a correct k-d tree after each insertion or erasure (default off)
 * 
 * -f Search for the tuple and the next tuple after erasure of each tuple (default off)
 * 
 * -w Create a worst-case set of coordinates by walking a k-d tree in order (default off)
 * 
 * -r Reverse the order of the worst-case set of coordinates for erasure (default off)
 *
 * -h Help
 */

 /*
  * Include kdTreeDynamic.h first so that KD_TREE_DYNAMIC_H will be defined
  * for kdTreeMergeSort.h, kdTreeKnlogn.h and kdTreeNlogn.h
  */
#include "kdTreeDynamic.h"

#ifdef NLOGN
#include "kdTreeNlogn.h"
#else
#include "kdTreeKnlogn.h"
#endif

/*
 * This is the type used for the test. Change the intrisic type in
 * this typedef to test the k-d tree with different intrisic types.
 */
typedef int64_t kdKey_t; // Add required #include and using to kdTreeNode.h

/*
  * Calculate the mean and standard deviation of the elements of a vector.
  *
  * Calling parameter:
  *
  * vec - a vector
  * 
  * return a pair that contains the mean and standard deviation
  */

 template <typename T>
 pair<double, double> calcMeanStd(vector<T> const& vec) {
  double sum = 0, sum2 = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    sum += vec[i];
    sum2 += vec[i] * vec[i];
  }
double n = static_cast<double>(vec.size());
return make_pair(sum / n, sqrt((n * sum2) - (sum * sum)) / n);
}

/* Create and search a k-d tree. */
int main(int argc, char** argv) {

  // Set the defaults then parse the input arguments.
  size_t iterations = 1;
  signed_size_t numPoints = 262144;
  signed_size_t numNeighbors = 5;
  signed_size_t numDimensions = 3;
  signed_size_t numThreads = 1;
  signed_size_t maximumNumberOfNodesToPrint = 5;
  kdKey_t searchDivisor = 10;
  bool balanced = false;
  bool neighbors = false;
  bool verify = false;
  bool find = false;
  bool worst = false;
  bool reverse = false;

  for (signed_size_t i = 1; i < argc; ++i) {
    if (0 == strcmp(argv[i], "-i") || 0 == strcmp(argv[i], "--iterations")) {
      iterations = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-n") || 0 == strcmp(argv[i], "--numPoints")) {
      numPoints = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-m") || 0 == strcmp(argv[i], "--numNeighbors")) {
      numNeighbors = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-d") || 0 == strcmp(argv[i], "--numDimensions")) {
      numDimensions = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-t") || 0 == strcmp(argv[i], "--numThreads")) {
      numThreads = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-s") || 0 == strcmp(argv[i], "--searchDivisor")) {
      searchDivisor = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-p") || 0 == strcmp(argv[i], "--maximumNodesToPrint")) {
      maximumNumberOfNodesToPrint = atol(argv[++i]);
      continue;
    }
    if (0 == strcmp(argv[i], "-b") || 0 == strcmp(argv[i], "--balanced")) {
      balanced = !balanced;
      continue;
    }
    if (0 == strcmp(argv[i], "-g") || 0 == strcmp(argv[i], "--neighbors")) {
      neighbors = !neighbors;
      continue;
    }
    if (0 == strcmp(argv[i], "-v") || 0 == strcmp(argv[i], "--verify")) {
      verify = !verify;
      continue;
    }
    if (0 == strcmp(argv[i], "-f") || 0 == strcmp(argv[i], "--find")) {
      find = !find;
      continue;
    }
    if (0 == strcmp(argv[i], "-w") || 0 == strcmp(argv[i], "--worst")) {
      worst = !worst;
      continue;
    }
    if (0 == strcmp(argv[i], "-r") || 0 == strcmp(argv[i], "--reverse")) {
      reverse = !reverse;
      continue;
    }

    if (0 == strcmp(argv[i], "-h") || 0 == strcmp(argv[i], "--help")) {
      cout << endl << "Usage:" << endl << endl
           << "kdTreeKnlogn [-n N] [-m M] [-x X] [-d D] [-t T] [-s S] [-p P] [-z Z] [-b] [-c] [-r] [-v] [-f]" << endl << endl
           << "where the command-line options are interpreted as follows." << endl << endl
           << "-i The number I of iterations of k-d tree creation" << endl << endl
           << "-n The number N of randomly generated points used to build the k-d tree" << endl << endl
           << "-m The maximum number M of nearest neighbors added to a priority queue" << endl << endl
           << "-d The number of dimensions D (aka k) of the k-d tree" << endl << endl
           << "-t The number of threads T used to build and search the k-d tree" << endl << endl
           << "-s The search divisor S used for region search" << endl << endl
           << "-p The maximum number P of nodes to report when reporting region search results" << endl << endl
           << "-b Build a balanced k-d tree for comparison to the dynamic k-d tree" << endl << endl
           << "-g Find nearest neighbors to each point" << endl << endl
           << "-v Verify the k-d tree ordering and balance after insertion or erasure of each point" << endl << endl
           << "-f Check for the next point after deleting each point (a cheap tree-order check)" << endl << endl
           << "-w Create a worst-case set of coordinates by walking a k-d tree in order" << endl << endl
           << "-r Reverse the order of the worst-case set of coordinates for erasure" << endl << endl
           << "-h List the command-line options" << endl << endl;
      exit(1);
    }
    {
      ostringstream buffer;
      buffer << "\n\nillegal command-line argument: " << argv[i] << endl;
      throw runtime_error(buffer.str());
    }
  }


  // Calculate the number of child threads to be the number of threads minus 1, then
  // calculate the maximum tree depth at which to launch a child thread.  Truncate
  // this depth such that the total number of threads, including the master thread, is
  // an integer power of 2, hence simplifying the launching of child threads by restricting
  // them to only the < branch of the tree for some depth in the tree.
  signed_size_t n = 0;
  if (numThreads > 0) {
    while (numThreads > 0) {
      ++n;
      numThreads >>= 1;
    }
    numThreads = 1 << (n - 1);
  }
  else {
    numThreads = 0;
  }
  signed_size_t const childThreads = numThreads - 1;
  signed_size_t maximumSubmitDepth = -1;
  if (numThreads < 2) {
    maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
  }
  else if (numThreads == 2) {
    maximumSubmitDepth = 0;
  }
  else {
    maximumSubmitDepth = static_cast<signed_size_t>(floor(log(static_cast<double>(childThreads)) / log(2.)));
  }
  cout << endl << "max number of threads = " << numThreads
       << "  max submit depth = " << maximumSubmitDepth << endl << endl;

  // Create an instance of KdTreeDynamic.
  auto tree = new KdTreeDynamic<kdKey_t>(maximumSubmitDepth);

  // Declare and initialize the coordinates and oneCoordinte vectors.
  vector<vector<kdKey_t>> coordinates(numPoints, vector<kdKey_t>(numDimensions));
  vector<kdKey_t> oneCoordinate(numPoints);

  // Calculate a delta coordinate by dividing the positive range of int64_t
  // by the number of points and truncating the quotient. Because the positive
  // range is less than half the full range of int64_t, multiplying the
  // delta coordinate by the number of points ought to produce a product
  // that is less than half the full range of int64_t and therefore avoid
  // possible overflow when comparing keys via the superKeyCompare function.
  // Calculate a padding coordinate to center the coordinates about zero.
  signed_size_t deltaCoordinate = LLONG_MAX / numPoints;
  size_t padCoordinate = (ULLONG_MAX - (numPoints * deltaCoordinate)) / 2;

  // Initialize each tuple. Equally space each coordinate
  // across the range of centered coordinates.
  kdKey_t beginCoordinate = LLONG_MIN + padCoordinate;
  kdKey_t thisCoordinate = beginCoordinate;
  kdKey_t endCoordinate = 0;
  for (signed_size_t i = 0; i < numPoints; ++i) {
    oneCoordinate[i] = thisCoordinate;
    endCoordinate = thisCoordinate;
    thisCoordinate += deltaCoordinate;
  }

  cout << "deltaCoordinate = " << deltaCoordinate << endl;
  cout << "padCoordinate = " << padCoordinate << endl;
  cout << "beginCoordinate = " << beginCoordinate << endl;
  cout << "endCoordinate = " << endCoordinate << endl << endl;

  // Allocate vectors to store the execution times and statistics.
  vector<double> createTime(iterations);
  vector<double> insertTime(iterations);
  vector<double> verifyTime(iterations);
  vector<double> eraseTime(iterations);
  vector<double> containsTime(iterations);
  vector<double> neighborsTimeStatic(iterations);
  vector<double> neighborsTimeDynamic(iterations);

  // If a worst-case set of coordinates is requested, create that set
  // by creating a static, balanced k-d tree and walking that tree in order.
  if (worst) {

    // Shuffle the coordinates vector independently for each dimension.
    std::mt19937_64 g(std::mt19937_64::default_seed);
    for (signed_size_t j = 0; j < numDimensions; ++j) {
      shuffle(oneCoordinate.begin(), oneCoordinate.end(), g);
      for (signed_size_t i = 0; i < numPoints; ++i) {
        coordinates[i][j] = oneCoordinate[i];
      }
    }

    // Create a static KdTree instance tree.
    signed_size_t numNodes;
    double allocateTime, sortTime, removeTime, kdTime,
            verifyTime, deallocateTime, unsortTime;
    auto const balancedTree =
      KdTree<kdKey_t>::createKdTree(coordinates, maximumSubmitDepth, numNodes,
                                    allocateTime, sortTime, removeTime, kdTime,
                                    verifyTime, deallocateTime, unsortTime);

    // Re-size the coordinates vector in case duplicate coordinates were removed.
    coordinates.resize(numNodes);

    // Create a KdTreeDynamic tree instance that has as its root node
    // the root node of the static KdTree instance. The KdTreeDynamic
    // constructor provides a way to create a static k-d tree that
    // can thereafter be modified as a dynamic k-d tree.
    auto dynamicTree = new KdTreeDynamic<kdKey_t>(maximumSubmitDepth,
                                                  balancedTree->getRoot());

    // Walk the static k-d tree in increasing order and
    // copy each tuple into a coordinate, which sorts the
    // coordinates in the coordinates vector.
    size_t count = dynamicTree->getSortedTree(coordinates);

    // Delete the static KdTree instance, which does not recursively
    // delete the KdNode instances when KD_TREE_DYNAMIC_H is defined
    // (see the ~KdTree destructor). Then delete the KdTreeDynamic
    // instance, which does delete the KdNode instances.
    delete balancedTree;
    delete dynamicTree;
  }

  // Iterate the construction of the k-d tree to improve statistics.
  signed_size_t numberOfNodes = 0, staticNumberOfNodes = 0;
  size_t treeHeight = 0, staticTreeHeight = 0;
  std::mt19937_64 g(std::mt19937_64::default_seed);
  for (size_t k = 0; k < iterations; ++k) {

    // Shuffle the coordinates vector independently for each dimension.
    // Skip this step if worst-case coordinates have already been created.
    if (!worst) {
      for (signed_size_t j = 0; j < numDimensions; ++j) {
        shuffle(oneCoordinate.begin(), oneCoordinate.end(), g);
        for (signed_size_t i = 0; i < numPoints; ++i) {
          coordinates[i][j] = oneCoordinate[i];
        }
      }
    }

    // Create a static, balanced k-d tree from the coordinates.
    if (balanced) {

      // Create the static k-d tree.
      signed_size_t numNodes;
      double allocateTime, sortTime, removeTime, kdTime,
             verifyTime, deallocateTime, unsortTime;
      KdTree<kdKey_t>* const tree =
        KdTree<kdKey_t>::createKdTree(coordinates, maximumSubmitDepth, numNodes,
                                      allocateTime, sortTime, removeTime, kdTime,
                                      verifyTime, deallocateTime, unsortTime);

      // Record the time for k-d tree creation, ignoring verifyTime and unsortTime.
      createTime[k] = allocateTime + sortTime + removeTime + kdTime + deallocateTime;

      // Record the number of nodes and the tree height for the static tree.
      staticNumberOfNodes = numNodes;
      staticTreeHeight = KdTreeDynamic<kdKey_t>::getHeight(tree->getRoot());

      // Find numNeighbors nearest neighbors to each coordinate.
      if (neighbors) {
        forward_list< pair<cpp_int, KdNode<kdKey_t>*> > neighborList;
        vector<kdKey_t> query(coordinates[0].size());
        auto beginTime = steady_clock::now();
        for (size_t i = 0; i < coordinates.size(); ++i) {
          neighborList.clear();
          for (size_t j = 0; j < query.size(); ++j) {
          query[j] = coordinates[i][j];
          }
        tree->findNearestNeighbors(neighborList, query, numNeighbors);
        }
        auto endTime = steady_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
        neighborsTimeStatic[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
      }

      // No further need for the balanced k-d tree, so delete it.
      delete tree;
    }

    // Insert each coordinate into the dynamic k-d tree.
    auto beginTime = steady_clock::now();
    for (size_t i = 0; i < coordinates.size(); ++i) {
      if (tree->insert(coordinates[i])) {
        if (verify) {
          tree->verifyKdTree(numDimensions, maximumSubmitDepth);
        }
      } else {
        cout << "\n\nfailed to insert tuple:";
        tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
        cout << endl;
        exit(0);
      }
    }
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    insertTime[k] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Verify correct order of each node in the k-d tree and count the nodes.
    beginTime = steady_clock::now();
    numberOfNodes = tree->verifyKdTree(numDimensions, maximumSubmitDepth);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    verifyTime[k] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
    if (static_cast<size_t>(numberOfNodes) != coordinates.size()) {
      ostringstream buffer;
      buffer << "\n\nnumber of coordinates = " << coordinates.size()
             << "  number of nodes = " << numberOfNodes << endl;
      throw runtime_error(buffer.str());
     } else {
      treeHeight = KdTreeDynamic<kdKey_t>::getHeight(tree->getRoot());
     }

    // Search for each coordinate in the k-d tree. 
    //
    // No need to reshuffle the coordinates prior to searching the tree
    // for each coordinate because search does not rebalance the tree;
    // hence, the insertion order of the coordinates is irrelevant to search.
    beginTime = steady_clock::now();
    for (size_t i = 0; i < coordinates.size(); ++i) {
      if (!tree->contains(coordinates[i])) {
        cout << "\n\nfailed to find tuple:";
        tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
        cout << endl;
        exit(0);
      }
    }
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    containsTime[k] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    // Find numNeighbors nearest neighbors to each coordinate.
    if (neighbors) {
      forward_list< pair<cpp_int, KdNode<kdKey_t>*> > neighborList;
      vector<kdKey_t> query(coordinates[0].size());
      auto beginTime = steady_clock::now();
      for (size_t i = 0; i < coordinates.size(); ++i) {
        neighborList.clear();
        for (size_t j = 0; j < query.size(); ++j) {
        query[j] = coordinates[i][j];
        }
      tree->findNearestNeighbors(neighborList, query, numNeighbors);
      }
      auto endTime = steady_clock::now();
      auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      neighborsTimeDynamic[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
    }

    // Reshuffle the coordinates prior to erasing each coordinate from
    // the tree because erasure rebalances the tree and hence the insertion
    // order of the keys may influence the performance of erasure. But
    // skip this step if worst-case coordinates have already been created.
    //
    // It is necessary to shuffle a vector of pointers to tuple arrays
    // because the shuffle function won't shuffle a 2D vector.
    if (!worst) {
      auto saveCoordinates = coordinates;
      auto tuplePointers = vector<kdKey_t*>(coordinates.size());
      for (size_t i = 0; i < coordinates.size(); ++i) {
        tuplePointers[i] = saveCoordinates[i].data();
      }
      shuffle(tuplePointers.begin(), tuplePointers.end(), g);
      for (size_t i = 0; i < coordinates.size(); ++i) {
        for (size_t j = 0; j < coordinates[0].size(); ++j) {
          coordinates[i][j] = tuplePointers[i][j];
        }
      }
    }

    // Erase each coordinate from the dynamic k-d tree,
    // and reverse the order of the coordinates if both
    // worst and reverse are true.
    beginTime = steady_clock::now();
    if (worst && reverse) {
      for (signed_size_t i = coordinates.size() - 1; i >= 0; --i) {
        if (tree->erase(coordinates[i])) {
          if (verify) {
            tree->verifyKdTree(numDimensions, maximumSubmitDepth);
          }
          if (find && tree->contains(coordinates[i])) {
            cout << "\n\nfound tuple after erasing tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
          if (find && i > 0 && !tree->contains(coordinates[i-1])) {
            cout << "\n\nfailed to find next tuple after erasing tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
        } else {
            cout << "\n\nfailed to erase tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
      }
    } else {
      for (size_t i = 0; i < coordinates.size(); ++i) {
        if (tree->erase(coordinates[i])) {
          if (verify) {
            tree->verifyKdTree(numDimensions, maximumSubmitDepth);
          }
          if (find && tree->contains(coordinates[i])) {
            cout << "\n\nfound tuple after erasing tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
          if (find && i < coordinates.size()-1 && !tree->contains(coordinates[i+1])) {
            cout << "\n\nfailed to find next tuple after erasing tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
        } else {
            cout << "\n\nfailed to erase tuple:";
            tree->printTuple(coordinates[i]); // Need to implement printTupleToStream
            cout << endl;
            exit(0);
        }
      }
    }
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    eraseTime[k] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    if ( !tree->isEmpty() ) {
      throw runtime_error("\n\ntree is not empty\n");
    }

    cout << "finished iteration " << (k + 1) << endl;
  }

  // Report the k-d tree statistics.
  cout << endl << "dynamic tree:" << endl << endl;
  cout << "number of nodes = " << numberOfNodes
               << "  k-d tree height = " << treeHeight << endl << endl;

  auto timePair = calcMeanStd<double>(insertTime);
  cout << "insert time = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd<double>(verifyTime);
  cout << "verify time = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd<double>(containsTime);
  cout << "search time = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd<double>(eraseTime);
  cout << "delete time = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  if (neighbors) {
    timePair = calcMeanStd<double>(neighborsTimeDynamic);
    cout << "neighbors time = " << fixed << setprecision(4) << timePair.first
        << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;
  }
  cout << endl;

  if (balanced) {
    cout << "static tree:" << endl << endl;
    cout << "number of nodes = " << staticNumberOfNodes
                << "  k-d tree height = " << staticTreeHeight << endl << endl;

    timePair = calcMeanStd<double>(createTime);
    cout << "create time = " << fixed << setprecision(4) << timePair.first
        << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;
    if (neighbors) {
      timePair = calcMeanStd<double>(neighborsTimeStatic);
      cout << "neighbors time = " << fixed << setprecision(4) << timePair.first
          << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;
    }
    cout << endl;
  }
  // Delete the k-d tree instance.
  delete tree;

  return 0;
}

