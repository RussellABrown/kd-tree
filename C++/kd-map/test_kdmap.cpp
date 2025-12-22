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
 * Test program for kdMapKnlogn.h and kdMapNlogn.h
 *
 * Compile via: g++ -O3 -std=c++20 -pthread -W -D PREALLOCATE test_kdmap.cpp
 *
 * Optional compilation defines are as follows.
 * 
 * -D NLOGN - Select the O(n log n) algorithm instead of the O(kn log n) algorithm.
 * 
 * The following compilation defines apply to both O(n log n) and O(kn log n) algorithms.
 * 
 * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
 *
 * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
 *                              in the MergeSort::mergeSort* functions (default 15)
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
 *                     the calculated median in KdNode::partition (default 512)
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
 * test_kdmap [-i I] [-n N] [-m M] [-x X] [-d D] [-t T] [-s S] [-p P] [-r]
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
 * -g Find nearest neighbors to a query point near the origin
 * 
 * -j Perform a region search in a hypercubed centered at a query point near the origin
 *
 * -r Construct nearest-neighbors and reverse-nearest-neighbors maps
 *
 * -h Help
 */

#include "kdMap.h"

#ifdef NLOGN
#include "kdMapNlogn.h"
#else
#include "kdMapKnlogn.h"
#endif

/*
 * These are the types used for the test. Change the intrisic types in
 * these typedefs to test the k-d tree with different intrisic types.
 */
typedef int64_t kdKey_t;  // Add required #include and using to kdMapNode.h
typedef string kdValue_t; // Add required #include and using to kdMapNode.h

/*
  * Calculate the mean and standard deviation of the elements of a vector.
  *
  * Calling parameter:
  *
  * vec - a vector
  * 
  * return a pair that contains the mean and standard deviation
  */
pair<double, double> calcMeanStd(const vector<double>& vec) {
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
  bool neighbors = false;
  bool region = false;
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
    if (0 == strcmp(argv[i], "-g") || 0 == strcmp(argv[i], "--neighbors")) {
      neighbors = !neighbors;
      continue;
    }
    if (0 == strcmp(argv[i], "-j") || 0 == strcmp(argv[i], "--region")) {
      region = !region;
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
           << "kdTreeKnlogn [-n N] [-m M] [-x X] [-d D] [-t T] [-s S] [-p P] [-g] [-j] [-r] [-h]" << endl << endl
           << "where the command-line options are interpreted as follows." << endl << endl
           << "-i The number I of iterations of k-d tree creation" << endl << endl
           << "-n The number N of randomly generated points used to build the k-d tree" << endl << endl
           << "-m The maximum number M of nearest neighbors added to a priority queue" << endl << endl
           << "-x The number X of duplicate points added to test removal of duplicate points" << endl << endl
           << "-d The number of dimensions D (aka k) of the k-d tree" << endl << endl
           << "-t The number of threads T used to build and search the k-d tree" << endl << endl
           << "-s The search divisor S used for region search" << endl << endl
           << "-p The maximum number P of nodes to report when reporting region search results" << endl << endl
           << "-b Build a balanced k-d tree for comparison to the dynamic k-d tree" << endl << endl
           << "-g Find nearest neighbors to a query point near the origin" << endl << endl
           << "-j Perform a region search in a hypercube centered at a query point near the origin" << endl << endl
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

  // It is impossible to find more nearest neighbors than there are points.
  numNeighbors = min(numNeighbors, numPoints + extraPoints + 1);

  // Declare and initialize the coordinates and oneCoordinte vectors. Each element
  // of the coordinates vector is a pair wherein first is a vector of (x, y, z, w...)
  // coordinates and second is the string representation of an integer.
  vector<kdKey_t> oneCoordinate(numPoints);
  extraPoints = (extraPoints < numPoints) ? extraPoints : numPoints - 1;

  vector<pair<vector<kdKey_t>, kdValue_t>> coordinates(numPoints + extraPoints);
  for (size_t i = 0; i < coordinates.size(); ++i) {
    ostringstream buffer;
    buffer << i;
    coordinates[i] = make_pair<vector<kdKey_t>, kdValue_t>(vector<kdKey_t>(numDimensions), buffer.str());
  }

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

  cout << endl << "deltaCoordinate = " << deltaCoordinate << endl;
  cout << "padCoordinate = " << padCoordinate << endl;
  cout << "beginCoordinate = " << beginCoordinate << endl;
  cout << "endCoordinate = " << endCoordinate << endl;

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
  vector<double> kdTotalTime(iterations);
  vector<double> containsTime(iterations);
  vector<double> neighborsTime(iterations);
  vector<double> bruteNeighborsTime(iterations);
  vector<double> regionTime(iterations);
  vector<double> bruteRegionTime(iterations);

  // Initialize the Mersenne twister pseudo-random number generator.
  std::mt19937_64 g(std::mt19937_64::default_seed);

  // Create query arrays for searching the k-d tree via region and nearest-neighbors searches.
  vector <kdKey_t>query(numDimensions);
  vector <kdKey_t>queryLower(numDimensions);
  vector <kdKey_t>queryUpper(numDimensions);
  for (signed_size_t i = 0; i < numDimensions; i++) {
      query[i] = i;
      queryLower[i] = query[i] + (beginCoordinate / searchDivisor);
      queryUpper[i] = query[i] + (endCoordinate / searchDivisor);
  }
  
  // Iterate the creation of the k-d tree to improve statistics.
  signed_size_t numberOfNodes = 0;
  size_t numRegionNodes = 0, numNeighborsNodes = 0;
  unique_ptr<KdTree<kdKey_t, kdValue_t>> tree;
  for (size_t k = 0; k < iterations; ++k) {

    // Shuffle the coordinates vector independently for each dimension.
    for (signed_size_t j = 0; j < numDimensions; ++j) {
      shuffle(oneCoordinate.begin(), oneCoordinate.end(), g);
      for (signed_size_t i = 0; i < numPoints; ++i) {
        coordinates[i].first[j] = oneCoordinate[i];
      }
    }

    // Reflect tuples across coordinates[numPoints - 1] to initialize the extra points.
    for (signed_size_t i = 1; i <= extraPoints; ++i) {
      for (signed_size_t j = 0; j < numDimensions; ++j) {
        coordinates[numPoints - 1 + i].first[j] = coordinates[numPoints - 1 - i].first[j];
      }
    }

    // Create a static k-d tree and manage it via the 'tree' unique_ptr.
    // NOTE: reset() deletes any static k-d tree previously managed by 'tree'.
    vector<pair<vector<kdKey_t>, kdValue_t>> copyCoordinates = coordinates;
    signed_size_t numNodes;
    double aT, sT, rT, kT, vT;
    tree.reset(KdTree<kdKey_t, kdValue_t>::createKdTree(copyCoordinates, maximumSubmitDepth,
                                                        numNodes, aT, sT, rT, kT, vT));

    // Record the time for k-d tree creation, ignoring verifyTime and unsortTime.
    kdTotalTime[k] += aT + sT + rT + kT;
    allocateTime[k] += aT;
    sortTime[k] += sT;
    removeTime[k] += rT;
    kdTime[k] += kT;
    verifyTime[k] += vT;
   
    // Record the number of nodes and the tree height for the static tree.
    numberOfNodes = numNodes;

    // Find numNeighbors nearest neighbors to each coordinate.
    //
    // NOTE: neighborList and bruteList are lists of KdNode*
    // not of shared_ptr<KdNode> so their KdNode* will remain valid
    // until the 'tree' unique_ptr goes out of scope.
    if (neighbors) {
      forward_list<pair<double, KdNode<kdKey_t, kdValue_t>*>> neighborList;
      auto beginTime = steady_clock::now();
      tree->findNearestNeighbors(neighborList, query, numNeighbors);
      auto endTime = steady_clock::now();
      auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      neighborsTime[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
      numNeighborsNodes = distance(neighborList.begin(), neighborList.end());

      forward_list<pair<double, KdNode<kdKey_t, kdValue_t>*>> bruteList;
      beginTime = steady_clock::now();
      tree->bruteNearestNeighbors(bruteList, query, numNeighbors);
      endTime = steady_clock::now();
      duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      bruteNeighborsTime[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

      // Compare the results of nearest-neighbor search and brute-force search.
      tree->verifyNearestNeighbors(neighborList, bruteList);
    }

    // Perform a region search within a hypercube centered near the origin.
    //
    // NOTE: fastRegionList and slowRegionList are lists of KdNode*
    // not of shared_ptr<KdNode> so their KdNode* will remain valid
    // until the 'tree' unique_ptr goes out of scope.
    if (region) {
      list<KdNode<kdKey_t, kdValue_t>*> fastRegionList;
      auto beginTime = steady_clock::now();
      tree->searchRegion(fastRegionList, queryLower, queryUpper, maximumSubmitDepth, true);
      auto endTime = steady_clock::now();
      auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      regionTime[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;
      numRegionNodes = fastRegionList.size();

      list<KdNode<kdKey_t, kdValue_t>*> slowRegionList;
      beginTime = steady_clock::now();
      tree->searchRegion(slowRegionList, queryLower, queryUpper, maximumSubmitDepth, false);
      endTime = steady_clock::now();
      duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
      bruteRegionTime[k] = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

      // Compare the results of region search and brute-force search.
      tree->verifyRegionSearch(fastRegionList, slowRegionList);
    }

    cout << "finished iteration " << (k + 1) << endl;
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

  timePair = calcMeanStd(kdTime);
  cout << "kdBuildTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  timePair = calcMeanStd(verifyTime);
  cout << "verifyTime = " << fixed << setprecision(4) << timePair.first
       << setprecision(4) << "  std dev = " << timePair.second << " seconds" << endl;

  if (neighbors) {
      timePair = calcMeanStd(neighborsTime);
      cout << "\nneighbors time = " << fixed << setprecision(6) << timePair.first
           << setprecision(6) << "  std dev = " << timePair.second << " seconds" << endl;
      timePair = calcMeanStd(bruteNeighborsTime);
      cout << "brute time     = " << fixed << setprecision(6) << timePair.first
          << setprecision(6) << "  std dev = " << timePair.second << " seconds" << endl;
  }

  if (region) {
      cout << "\nFound " << numRegionNodes << " tuples by region search within "
          << (queryUpper[0] - queryLower[0]) << " units of ";
      tree->printTuple(query);
      cout << " in all dimensions.\n" << endl;
      timePair = calcMeanStd(regionTime);
      cout << "region time = " << fixed << setprecision(6) << timePair.first
          << setprecision(6) << "  std dev = " << timePair.second << " seconds" << endl;
      timePair = calcMeanStd(bruteRegionTime);
      cout << "brute time  = " << fixed << setprecision(6) << timePair.first
          << setprecision(6) << "  std dev = " << timePair.second << " seconds" << endl;
 }

  cout << endl;

#ifdef REVERSE_NEAREST_NEIGHBORS
  // Optionally construct a nearest neighbor vector and a reverse nearest neighbors vector.
  // Each vector element contains a list that is initialized to an empty list.
  //
  // NOTE: the nearest neighbor vector and reverse nearest neighbor vector contain
  // lists of KdNode* not of shared_ptr<KdNode> so their KdNode* will remain value
  // until the 'tree' unique_ptr goes out of scope.
  if (reverseNearestNeighbors) {
    auto beginTime = steady_clock::now();
    vector<forward_list<pair<double, KdNode<kdKey_t, kdValue_t>*>>> nn(coordinates.size());
    vector<forward_list<pair<double, KdNode<kdKey_t, kdValue_t>*>>> rnn(coordinates.size());
    vector<mutex> mutexes(coordinates.size());
    auto endTime = steady_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const vectorTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "vector initialization  time = " << fixed << setprecision(3) << vectorTime << " seconds" << endl << endl;

    beginTime = steady_clock::now();
    tree->findReverseNearestNeighbors(nn, rnn, mutexes, numDimensions, numNeighbors, maximumSubmitDepth);
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const reverseNearestNeighborTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "reverse nearest neighbor time = " << fixed << setprecision(3) << reverseNearestNeighborTime << " seconds" << endl << endl;
    cout << "number of non-empty nearest neighbors lists = " << tree->nonEmptyLists(nn) << endl;
    cout << "number of non-empty reverse nearest neighbors lists = " << tree->nonEmptyLists(rnn) << endl << endl;

    // Report the mean and standard deviation distance and number of reverse nearest neighbors.
    double meanSize, stdSize, meanDist, stdDist;
    tree->calculateMeanStd(rnn, meanSize, stdSize, meanDist, stdDist);
    cout << "mean reverse nearest neighbor distance = " << scientific << setprecision(3) << meanDist
         << "  standard deviation = " << stdDist << endl;
    cout << "mean reverse nearest neighbor list size = " << fixed << setprecision(3) << meanSize
         << "  standard deviation = " << stdSize << endl << endl;

    // Report the mean and standard deviation distance and number of nearest neighbors.
    tree->calculateMeanStd(nn, meanSize, stdSize, meanDist, stdDist);
    cout << "mean nearest neighbor distance = " << scientific << meanDist
         << "  standard deviation = " << stdDist << endl;
    cout << "mean nearest neighbor list size = " << fixed << setprecision(3) << meanSize
         << "  standard deviation = " << stdSize << endl << endl;

    // Verify the consistency between the nearest neighbors and reverse nearest neighbors vectors.
    beginTime = steady_clock::now();
    tree->verifyReverseNeighbors(nn, rnn, coordinates.size());
    endTime = steady_clock::now();
    duration = duration_cast<std::chrono::microseconds>(endTime - beginTime);
    double const verifyReverseTime = static_cast<double>(duration.count()) / MICROSECONDS_TO_SECONDS;

    cout << "verify reverse nearest neighbor time = " << fixed << setprecision(2) << verifyReverseTime << " seconds" << endl << endl;
  }
#endif

  return 0;
}
