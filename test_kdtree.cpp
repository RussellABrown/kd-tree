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
 * Test program for kdTreeKnlogn.h, kdTreeNlogn.h, and kdTreeYuCao.h
 *
 * Compile via: g++ -O3 -std=c++11 -pthread -W -D PREALLOCATE test_kdtree.cpp
 *
 * Optional compilation defines are as follows.
 * 
 * -D NLOGN - Select the O(n log n) algorithm instead of the O(kn log n) algorithm.
 * 
 * -D YUCAO - Select the O(kn log n) presort and O(n log n) k-d tree build algorithms
 *            if -D NLOGN is not also specified. See the journal article by Yu Cao et al.,
 *            "A New Method to Construct the KD Tree Based on Presorted Results", Complexity
 *             Volume 2020, Article ID 8883945, https://doi.org/10.1155/2020/8883945
 * 
 * -D DEBUG_PRINT - If -D YUCAO is also specified, provide a simple coordinates vector
 *                  as specified in Figure 1 of the Yu Cao et al. article and test it
 *                  via this test_kdTree.cpp file and print various arrays the createKdTree
 *                  function of kdTreeYuCao.h to facilitate debugging Yu Cao's algorithm. * 
 * 
 * The following compilation defines create a k-d tree with more nodes than shown
 * in Figure 1 of the Yu Cao et al. article.
 * 
 * -D YUCAO -D DEBUG_PRINT -D LARGE_TREE
 * 
 * The following compilation defines create a k-d tree with an even number of nodes.
 * 
 * -D YUCAO -D DEBUG_PRINT-D LARGE_TREE -D EVEN_TREE
 * 
 * The following compilation defines apply to the O(n log n), O(kn log n) and Yu Cao algorithms.
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
 *  The following compilation define applies only to the O(kn log n) algorithm.
 * 
 * -D PARTITION_CUTOFF=n - A cutoff for switching from 1 to 2 threads to partition
 *                         reference arrays in KdTree::buildKdTree (default 4096)
 * 
 * The following compilation define applies only to the Yu Cao algorithm.
 * 
 * -D DUAL_THREAD_YUCAO - Execute the dual-threaded version of the k-d tree-building
 *                        algorithm depicted in Figure 2 of Yu Cao et al.
 * 
 * Usage:
 *
 * test_kdtree [-i I] [-n N] [-m M] [-x X] [-d D] [-t T] [-s S] [-p P] [-b] [-c] [-r]
 *
 * where the command-line options are interpreted as follows.
 * 
 * -i The number I of iterations of k-d tree creation
 *
 * -n The number N of randomly generated points used to build the k-d tree
 *
 * -m The maximum number M of nearest neighbors added to a priority queue
 *    when searching the k-d tree for nearest neighbors (default 5)
 *
 * -x The number X of duplicate points added to test removal of duplicate points
 *
 * -d The number of dimensions D (aka k) of the k-d tree
 *
 * -t The number of threads T used to build and search the k-d tree
 *
 * -s The search divisor S used for region search (default 10)
 *
 * -p The maximum number P of nodes to report when reporting region search results
 *
 * -b Compare k-d tree nearest neighbors search to exhaustive search
 *
 * -c Compare k-d tree region search to exhaustive search
 *
 * -r Construct nearest-neighbors and reverse-nearest-neighbors maps
 *
 * -h Help
 */

#ifdef NLOGN
#include "kdTreeNlogn.h"
#else
#ifdef YUCAO
#include "kdTreeYuCao.h"
#else
#include "kdTreeKnlogn.h"
#endif
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
pair<double, double> calcMeanStd(vector<double> const& vec) {
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
  signed_size_t extraPoints = 100;
  signed_size_t numDimensions = 3;
  signed_size_t numThreads = 4;
  signed_size_t maximumNumberOfNodesToPrint = 5;
  kdKey_t searchDivisor = 10;
  bool bruteForceSearch = false;
  bool bruteForceRegion = false;
  bool reverseNearestNeighbors = false;

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
    if (0 == strcmp(argv[i], "-x") || 0 == strcmp(argv[i], "--extraPoints")) {
      extraPoints = atol(argv[++i]);
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
    if (0 == strcmp(argv[i], "-b") || 0 == strcmp(argv[i], "--bruteForceSearch")) {
      bruteForceSearch = !bruteForceSearch;
      continue;
    }
    if (0 == strcmp(argv[i], "-c") || 0 == strcmp(argv[i], "--bruteForceRegion")) {
      bruteForceRegion = !bruteForceRegion;
      continue;
    }
#ifdef REVERSE_NEAREST_NEIGHBORS
    if (0 == strcmp(argv[i], "-r") || 0 == strcmp(argv[i], "--reverseNearestNeighbors")) {
      reverseNearestNeighbors = !reverseNearestNeighbors;
      continue;
    }
#endif
    if (0 == strcmp(argv[i], "-h") || 0 == strcmp(argv[i], "--help")) {
      cout << endl << "Usage:" << endl << endl
           << "kdTreeKnlogn [-n N] [-m M] [-x X] [-d D] [-t T] [-s S] [-p P] [-b] [-c] [-r] [-h]" << endl << endl
           << "where the command-line options are interpreted as follows." << endl << endl
           << "-i The number I of iterations of k-d tree creation" << endl << endl
           << "-n The number N of randomly generated points used to build the k-d tree" << endl << endl
           << "-m The maximum number M of nearest neighbors added to a priority queue" << endl << endl
#ifndef YUCAO
           << "-x The number X of duplicate points added to test removal of duplicate points" << endl << endl
#endif
           << "-d The number of dimensions D (aka k) of the k-d tree" << endl << endl
           << "-t The number of threads T used to build and search the k-d tree" << endl << endl
           << "-s The search divisor S used for region search" << endl << endl
           << "-p The maximum number P of nodes to report when reporting region search results" << endl << endl
           << "-b Compare k-d tree nearest neighbors search to exhaustive search" << endl << endl
           << "-c Compare k-d tree region search to exhaustive search" << endl << endl
#ifdef REVERSE_NEAREST_NEIGHBORS
           << "-r Construct nearest-neighbors and reverse-nearest-neighbors maps" << endl << endl
#endif
           << "-h List the command-line options" << endl << endl;
      exit(1);
    }
    {
      ostringstream buffer;
      buffer << "\n\nillegal command-line argument: " << argv[i] << endl;
      throw runtime_error(buffer.str());
    }
  }

  // Declare and initialize the coordinates and oneCoordinte vectors.
  extraPoints = (extraPoints < numPoints) ? extraPoints : numPoints - 1;
  vector<kdKey_t> oneCoordinate(numPoints);
#if !defined (YUCAO) || !defined(DEBUG_PRINT)
  vector<vector<kdKey_t>> coordinates(numPoints + extraPoints, vector<kdKey_t>(numDimensions));
#else
  // Override the coordinates vector with a definition from Yu Cao Figure 1.
#ifdef LARGE_TREE
#ifndef EVEN_TREE
  vector<vector<kdKey_t>> coordinates = { {2,3,4}, {5,4,2}, {9,6,7}, {4,7,9}, {8,1,5}, 
                                          {7,2,6}, {9,4,1}, {8,3,2}, {9,7,8}, {6,3,2},
                                          {3,4,5}, {1,6,8}, {9,5,3}, {2,1,3}, {8,7,5} };
#else
  vector<vector<kdKey_t>> coordinates = { {2,3,4}, {5,4,2}, {9,6,7}, {4,7,9}, {8,1,5}, 
                                          {7,2,6}, {9,4,1}, {8,3,2}, {9,7,8}, {6,3,2},
                                          {3,4,5}, {1,6,8}, {9,5,3}, {2,1,3} };
#endif
#else
  vector<vector<kdKey_t>> coordinates = { {6,1}, {9,4}, {1,2}, {3,5}, {2,8}, {7,9}, {8,3} };
#endif
  numDimensions = coordinates[0].size();
#endif

  // Calculate a delta coordinate by dividing the positive range of int64_t
  // by the number of points and truncating the quotient. Because the positive
  // range is less than half the full range of int64_t, multiplying the
  // delta coordinate by the number of points ought to produce a product
  // that is less than half the full range of int64_t and therefore avoid
  // possible overflow when comparing keys via the superKeyCompare function.
  // Calculate a padding coordinate to center the coordinates about zero.
  signed_size_t deltaCoordinate = LLONG_MAX / numPoints;
  size_t padCoordinate = (ULLONG_MAX - (numPoints * deltaCoordinate)) / 2;

  // Initialize each tuple excluding the extra points. Equally space each coordinate
  // across the range of centered coordinates.
  kdKey_t beginCoordinate = LLONG_MIN + padCoordinate;
  kdKey_t thisCoordinate = beginCoordinate;
  kdKey_t endCoordinate = 0;
  for (signed_size_t i = 0; i < numPoints; ++i) {
    oneCoordinate[i] = thisCoordinate;
    endCoordinate = thisCoordinate;
    thisCoordinate += deltaCoordinate;
  }

  // These two coordinates indicate whether integer overlow occurs for region search.
  size_t maxCoordinate = static_cast<size_t>(endCoordinate) * 2;
  size_t minCoordinate = static_cast<size_t>(-beginCoordinate) * 2;

  cout << endl << "deltaCoordinate = " << deltaCoordinate << endl;
  cout << "padCoordinate = " << padCoordinate << endl;
  cout << "beginCoordinate = " << beginCoordinate << endl;
  cout << "endCoordinate = " << endCoordinate << endl;
  cout << "minCoordinate = -" << minCoordinate << endl;
  cout << "maxCoordinate = " << maxCoordinate << endl;

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
  cout << endl << "Max number of threads = " << numThreads << "  max submit depth = "
       << maximumSubmitDepth << endl << endl;

  // Allocate vectors to store the execution times for k-d tree creation.
  vector<double> allocateTime(iterations);
  vector<double> sortTime(iterations);
  vector<double> removeTime(iterations);
  vector<double> kdTime(iterations);
  vector<double> verifyTime(iterations);
  vector<double> deallocateTime(iterations);
  vector<double> unsortTime(iterations);
  vector<double> kdTotalTime(iterations);

  // Iterate the creation of the k-d tree to improve statistics.
  signed_size_t numberOfNodes = 0;
  KdTree<kdKey_t>* root = nullptr;
  std::mt19937_64 g(std::mt19937_64::default_seed);
  for (size_t k = 0; k < iterations; ++k) {

    // Shuffle the coordinates vector independently for each dimension.
#if !defined (YUCAO) || !defined(DEBUG_PRINT)
    for (signed_size_t j = 0; j < numDimensions; ++j) {
      shuffle(oneCoordinate.begin(), oneCoordinate.end(), g);
      for (signed_size_t i = 0; i < numPoints; ++i) {
        coordinates[i][j] = oneCoordinate[i];
      }
    }

    // Reflect tuples across coordinates[numPoints - 1] to initialize the extra points.
    for (signed_size_t i = 1; i <= extraPoints; ++i) {
      for (signed_size_t j = 0; j < numDimensions; ++j) {
        coordinates[numPoints - 1 + i][j] = coordinates[numPoints - 1 - i][j];
      }
    }
#endif

    // Create the k-d tree and record the execution times.
    root = KdTree<kdKey_t>::createKdTree(coordinates, maximumSubmitDepth, numberOfNodes,
                                         allocateTime[k], sortTime[k], removeTime[k], kdTime[k],
                                         verifyTime[k], deallocateTime[k], unsortTime[k]);
    kdTotalTime[k] =
      allocateTime[k] + sortTime[k] + removeTime[k] + kdTime[k] + verifyTime[k] + deallocateTime[k] + unsortTime[k];

#if !defined (YUCAO) || !defined(DEBUG_PRINT)
    // Verify that the k-d tree contains the correct number of nodes.
    if (numberOfNodes != numPoints) {
      ostringstream buffer;
      buffer << "\n\nk-d tree size = " << numberOfNodes
             << "  !=  number of points = " << numPoints << endl;
      throw runtime_error(buffer.str());
    }
#endif

    // If this iteration is not the final iteration, delete the root that
    // recursively deletes the k-d tree, including all instances of KdNode
    // and tuples arrays that were not deleted by KdNode::removeDuplicates.
    if (k < iterations - 1) {
      auto beginTime = steady_clock::now();
      delete root;
      auto endTime = steady_clock::now();
      auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      deallocateTime[k] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
    }
    cout << "finished iteration " << (k + 1) << endl;
  }

  // Allocate vectors to store execution times for k-d tree search.
  size_t searchIterations = min(iterations, static_cast<size_t>(numberOfNodes));
  vector<double> regionTime(searchIterations);
  vector<double> neighborTime(searchIterations);

  // Search the k-d tree via region search for the KdNodes that lie within a hyper-cube
  // that is centered at the first searchIterations (x, y, z, w...) coordinates and that
  // is (endCoordinate - beginCoordinate) wide in each dimension.

  vector<kdKey_t> query(numDimensions);
  vector<kdKey_t> queryLower(numDimensions);
  vector<kdKey_t> queryUpper(numDimensions);
  list<KdNode<kdKey_t>*> regionList;
  for (size_t i = 0; i < searchIterations; ++i) {
    for (signed_size_t j = 0; j < numDimensions; ++j) {
        query[j] = coordinates[i][j];
        queryLower[j] = query[j] + (beginCoordinate / searchDivisor);
        queryUpper[j] = query[j] + (endCoordinate / searchDivisor);
    }
    regionList.clear();
    auto beginTime = steady_clock::now();
    root->searchRegion(regionList, queryLower, queryUpper, maximumSubmitDepth);
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    regionTime[i] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
  }

  // It is impossible to find more nearest neighbors than there are points.
  numNeighbors = min(numNeighbors, numPoints + extraPoints + 1);

  // Search the k-d tree for up to numNeighbors nearest neighbors to the
  // first searchIterations (x, y, z, w...) coordinates.
  forward_list< pair<double, KdNode<kdKey_t>*> > neighborList;
  for (size_t i = 0; i < searchIterations; ++i) {
    for (signed_size_t j = 0; j < numDimensions; ++j) {
        query[j] = coordinates[i][j];
    }
    neighborList.clear();
    auto beginTime = steady_clock::now();
    root->findNearestNeighbors(neighborList, query, numNeighbors);
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    neighborTime[i] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
  }

  // Report the k-d tree statistics including creation and search execution times.
  cout << endl << "Number of nodes = " << numberOfNodes << endl << endl;

  auto timePair = calcMeanStd(kdTotalTime);
  cout << "kdTotalTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(allocateTime);
  cout << "allocateTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(sortTime);
  cout << "sortTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(removeTime);
  cout << "removeTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

#if defined(YUCAO) && !defined(NLOGN)
  timePair = calcMeanStd(unsortTime);
  cout << "unsortTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;
#endif

  timePair = calcMeanStd(kdTime);
  cout << "kdBuildTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;
  timePair = calcMeanStd(verifyTime);
  cout << "verifyTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(deallocateTime);
  cout << "deallocateTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl << endl;

  timePair = calcMeanStd(regionTime);
  cout << "regionTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(6) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(neighborTime);
  cout << "neighborTime = " << fixed << setprecision(6) << timePair.first
       << setprecision(7) << "  std dev = " << timePair.second << " seconds" << endl << endl;

  // Search the k-d tree via region search and brute force for the KdNodes
  // that lie within a hyper-cube centered near the origin.
  if (bruteForceRegion) {
    // Search the k-d tree via region search for the KdNodes that lie within a hyper-cube that is centered
    // near the origin and that extends +/- the search multiple times the delta coordinate.
    for (signed_size_t i = 0; i < numDimensions; ++i) {
      query[i] = i;
      queryLower[i] = query[i] + (beginCoordinate / searchDivisor);
      queryUpper[i] = query[i] + (endCoordinate / searchDivisor);
    }
    kdKey_t queryRange = (endCoordinate - beginCoordinate) / searchDivisor / 2;
    list<KdNode<kdKey_t>*> regionFast;
    root->searchRegion(regionFast, queryLower, queryUpper, maximumSubmitDepth);

    cout << regionFast.size() << " nodes within " << queryRange << " units of ";
    root->printTuple(query);
    cout << " in all dimensions." << endl << endl;
    cout << "List of the nearest <= " << maximumNumberOfNodesToPrint << " tree region-search k-d nodes within a "
        << queryRange << "-unit search distance follows:" << endl << endl;
    auto printRegionFast = root->sortByDistance(regionFast, query, maximumNumberOfNodesToPrint);
    root->printTuples(printRegionFast, maximumNumberOfNodesToPrint, numDimensions);
    cout << endl;

    // Verify that no duplicate KdNodes exist on the list returned from region search.
    auto itr1 = regionFast.begin();
    auto itr2 = itr1;
    ++itr2;
    for ( ; itr2 != regionFast.end(); ++itr1, ++itr2) {
      if (*itr1 == *itr2) {
        throw runtime_error("\n\nduplicate KdNode* on region-search list\n");
      }
    }

    auto beginTime = steady_clock::now();
    list<KdNode<kdKey_t>*> regionSlow;
    root->bruteRegion(regionSlow, queryLower, queryUpper);
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const slowRegionTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "brute-force region time = " << fixed << setprecision(6) << slowRegionTime << " seconds" << endl << endl;

    cout << regionSlow.size() << " nodes within " << queryRange << " units of ";
    root->printTuple(query);
    cout << " in all dimensions." << endl << endl;
    cout << "List of the nearest <= " << maximumNumberOfNodesToPrint << " brute region-search k-d nodes within a "
         << queryRange << "-unit search distance follows:" << endl << endl;
    auto printRegionSlow = root->sortByDistance(regionFast, query, maximumNumberOfNodesToPrint);
    root->printTuples(printRegionSlow, maximumNumberOfNodesToPrint, numDimensions);
    cout << endl;

    // Print the tree and brute distances.
    cout << "tree and brute closest region-search distances follow in increasing order:" << endl << endl;
    auto itf = printRegionFast.begin();
    auto its = printRegionSlow.begin();
    for ( ; itf != printRegionFast.end(); ++itf, ++its) {
      cout << scientific << sqrt(itf->first) << "\t" << sqrt(its->first) << endl;
    }
    cout << endl;

    // Verify that the region-search and brute-force lists are identical.
    // Both lists must be sorted before the KdNode* comparisons are
    // performed below because the region search and brute-force search
    // algorithms do not prepend KdNode* to their lists in the same order.
    regionFast.sort();
    regionSlow.sort();
    auto itrf = regionFast.begin();
    for (auto itrs = regionSlow.begin(); itrs != regionSlow.end(); ++itrf, ++itrs) {
      if (*itrf != *itrs) {
        throw runtime_error("\n\nnon-identical region-search and brute-force lists\n");
      }
    }
  }

  // Find nearest neighbors via nearest-neighbor search and brute force if requested.
  if (bruteForceSearch) {
    // Search the k-d tree for up to numNeighbors nearest neighbors to the query tuple.
    for (signed_size_t i = 0; i < numDimensions; ++i) {
      query[i] = i;
    }
    forward_list< pair<double, KdNode<kdKey_t>*> > neighborsFast;
    root->findNearestNeighbors(neighborsFast, query, numNeighbors);

    cout << "tree nearest-neighbor list size = " << distance(neighborsFast.begin(), neighborsFast.end()) << endl << endl;

    cout << "List of the nearest <= " << maximumNumberOfNodesToPrint << " tree neighbor-search k-d nodes follows:" << endl << endl;
    root->printTuples(neighborsFast, maximumNumberOfNodesToPrint, numDimensions);
    cout << endl;

    auto beginTime = steady_clock::now();
    forward_list< pair<double, KdNode<kdKey_t>*> > neighborsSlow;
    // Find only the number of nearest neighbors returned by findNearestNeighbors above.
    root->bruteNearestNeighbors(neighborsSlow, query, distance(neighborsFast.begin(), neighborsFast.end()));
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const slowNeighborTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "brute nearest-neighbor time = " << fixed << setprecision(6) << slowNeighborTime << " seconds" << endl << endl;
    cout << "brute neighbor list size = " << distance(neighborsSlow.begin(), neighborsSlow.end()) << endl << endl;

    cout << "List of the nearest <= " << maximumNumberOfNodesToPrint << " brute neighbor-search k-d nodes follows:" << endl << endl;
    root->printTuples(neighborsSlow, maximumNumberOfNodesToPrint, numDimensions);
    cout << endl;

    // Verify the consistency between the nearest neighbors lists
    // found by k-d tree search and by brute force.
    root->verifyNearestNeighbors(neighborsFast, neighborsSlow);

    // Print the tree and brute distances.
    cout << "tree and brute closest nearest-neighbor distances follow in increasing order:" << endl << endl;
    auto itf = neighborsFast.begin();
    auto its = neighborsSlow.begin();
    for ( ; itf != neighborsFast.end(); ++itf, ++its) {
      cout << scientific << sqrt(itf->first) << "\t" << sqrt(its->first) << endl;
    }
    cout << endl;
  }

#ifdef REVERSE_NEAREST_NEIGHBORS
  // Optionally construct a nearest neighbor vector and a reverse nearest neighbors vector.
  // Each vector element contains a list that is initialized to an empty list.
  if (reverseNearestNeighbors) {
    auto beginTime = steady_clock::now();
    vector< forward_list< pair<double, KdNode<kdKey_t>*> > > nn(coordinates.size());
    vector< forward_list< pair<double, KdNode<kdKey_t>*> > > rnn(coordinates.size());
    vector<mutex> mutexes(coordinates.size());
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const vectorTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "vector initialization  time = " << fixed << setprecision(2) << vectorTime << " seconds" << endl << endl;

    beginTime = steady_clock::now();
    root->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions, numNeighbors, maximumSubmitDepth);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const reverseNearestNeighborTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "reverse nearest neighbor time = " << fixed << setprecision(2) << reverseNearestNeighborTime << " seconds" << endl << endl;
    cout << "number of non-empty nearest neighbors lists = " << root->nonEmptyLists(nn) << endl;
    cout << "number of non-empty reverse nearest neighbors lists = " << root->nonEmptyLists(rnn) << endl << endl;

    // Report the mean and standard deviation distance and number of reverse nearest neighbors.
    double meanSize, stdSize, meanDist, stdDist;
    root->calculateMeanStd(rnn, meanSize, stdSize, meanDist, stdDist);
    cout << "mean reverse nearest neighbor distance = " << scientific << meanDist
         << "  standard deviation = " << stdDist << endl;
    cout << "mean reverse nearest neighbor list size = " << fixed << setprecision(3) << meanSize
         << "  standard deviation = " << stdSize << endl << endl;

    // Report the mean and standard deviation distance and number of nearest neighbors.
    root->calculateMeanStd(nn, meanSize, stdSize, meanDist, stdDist);
    cout << "mean nearest neighbor distance = " << scientific << meanDist
         << "  standard deviation = " << stdDist << endl;
    cout << "mean nearest neighbor list size = " << fixed << setprecision(3) << meanSize
         << "  standard deviation = " << stdSize << endl << endl;

    // Verify the consistency between the nearest neighbors and reverse nearest neighbors vectors.
    beginTime = steady_clock::now();
    root->verifyReverseNeighbors(nn, rnn, numberOfNodes);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const verifyReverseTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "verify reverse nearest neighbor time = " << fixed << setprecision(2) << verifyReverseTime << " seconds" << endl << endl;
  }
#endif

  // Delete the root that recursively deletes the k-d tree, including all
  // instances of KdNode and tuples arrays that were not deleted by the
  // KdNode::removeDuplicates function.
  auto beginTime = steady_clock::now();
  delete root;
  auto endTime = steady_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
  deallocateTime[iterations - 1] += static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

  return 0;
}
