/*
 * Copyright (c) 2015, 2019, 2020, 2023, 2025 Russell A. Brown
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
 * Test program for KdTree.java, KdTreeNlogn.java, KdTreeKnlogn.java, KdNode.java
 * MergeSort.java, NearestNeighborList.java, Pair.java, and Constants.java
 *
 * Configuration is controlled via the following constants in Constants.java
 *
 * NLOGN - If true, specifies the O[n log(n)] algorithm instead of the O[kn log(n)] algorithm.
 * 
 * ENABLE_PREFERRED_TEST - If true, compare the heights of a deleted 2-child node's
 *                          child subtrees to select a preferred replacement node.
 * 
 * ENABLE_1TO3 - If true, curtails recursive deletion when a subtree contains <= 3 nodes.
 * 
 * AVL_BALANCE - If true, the KdTreeDynamic.isBalanced method checks AVL balance;
 *               otherwise, this method checks for red-black balance.
 * 
 * HEIGHT_DIFF - An integer that specifies, for red-black balance, the maximum allowed
 *               height difference between the < and > subtrees of a node when one subtree
 *               is empty; or for AVL balance, the maximum allowed height difference
 *               between the < and > subtrees of a node. This balance is used to maintain
 *               the balance of a dynamic k-d tree, not a static k-d tree.
 * 
 * INSERTION_SORT_CUTOFF - An integer that specifies the number of nodes at which insertion
 *                         sort replaces merge sort in the MergeSort.mergeSort* methods.
 * 
 * MEDIAN_OF_MEDIANS_CUTOFF - An integer that specifies number of nodes at which insertion
 *                            sort replaces the median-of-medians algorithm in the
 *                            KdTreeNlogn.partiiton method.
 * 
 * KD_MAP_DYNAMIC - If true, the height of a node, which is defined as the maximum
 *                  number of nodes from that node to the bottom of the tree, is
 *                  computed. This height is used in the construction of a dynamic
 *                  k-d tree, but not in the construction of a static k-d tree.
 * 
 * MULTI_THREAD_CUTOFF - An integer that specifies the number of nodes above which
 *                       multiple threads are used to construct a static k-d tree.
 * 
 * NLOGN_CUTOFF - An integer the specifies the number of nodes above which multiple
 *                threads are used to construct a static k-d tree via the O[n log(n)]
 *                algorithm, independent of the number of threads specified by the
 *                -t command-line option (see below).
 * 
 * KNLOGN_CUTOFF - An integer the specifies the number of nodes above which multiple
 *                threads are used to construct a static k-d tree via the O[kn log(n)]
 *                algorithm, independent of the number of threads specified by the
 *                -t command-line option (see below).
 * 
 * MERGE_CUTOFF - An integer the specifies the number of nodes above which multiple
 *               threads are used for merge sort, independent of the number of threads
 *               specified by the -t command-line option (see below).
 * 
 * ENABLE_LINKED_LIST - If true, use the LinkedList version of region search
 *                      instead of the ArrayList version.
 * 
 * 
 * Usage:
 *
 * "java TestKdTreeDynamic [-n N] [-x X] [-d D] [-t T] [-c C] [-b] [-g] [-m M] \
 *                         [-j] [-s S] [-p P] [-i] [-h]
 *
 * where the command-line options are interpreted as follows.
 * 
 * -n The number N of randomly generated points used to build the k-d tree (default 262144)
 *
 * -x The number X of duplicate points added to the set of randomly generated points (default 100)
 *
 * -d The number of dimensions D (aka k) of the k-d tree (default 3)
 *
 * -t The number of threads T used to build and search the k-d tree (default 1)
 *
 * -c The multi-thread cutoff below which below which only a single thread is used
 *    to build the tree (default 65536)
 *
 * -b Build a balanced k-d tree for comparison to the dynamic k-d tree (default off)
 * 
 * -g Find nearest neighbors to a query point
 *
 * -m The maximum number M of nearest neighbors added to a priority queue
 *    when searching the k-d tree for nearest neighbors (default 5)
 * 
 * -j Perform a region search in a hypercube centered at a query point
 * 
 * -s The search divisor S that modifies the size of the hypercube for region search (default 10)
 *
 * -p The maximum number P of nodes to report when reporting nearest neighbor results (default 5)
 * 
 * -i The number I of iterations of k-d tree insertion, search, and deletion (default 1)
 *
 * -h Help, i.e., print the above explanation of command-line options.
 */
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.List;
import java.util.LinkedList;
import java.util.Random;

/**
 * <p>
 * The k-d tree is described by Jon Bentley in "Multidimensional Binary Search Trees Used
 * for Associative Searching", CACM 18(9): 509-517, 1975.  One approach to building a
 * balanced k-d tree involves finding the median of the data at each level of recursive
 * subdivision of those data.  A O(n) algorithm that may be used to find the median and
 * partition data about that median is described by Manuel Blum, Robert W. Floyd, Vaughan
 * Pratt, Ronald L. Rivest and Robert E. Tarjan in "Time Bounds for Selection", Journal of
 * Computer and System Sciences 7: 448-461, 1973 as well as by Thomas H. Cormen, Charles E.
 * Leiserson, Ronald L. Rivest and Clifford Stein in Section 9.3, "Selection in worst-case
 * linear time", pages 220-222 in "Introduction to Algorithms, Third Edition", MIT Press,
 * Cambridge, Massachusetts, 2009.
 * </p>
 *
 * @author Russell A. Brown
 */
public class TestKdTree {

    /**
     * <p>
     * The {@code shuffleArray} method shuffles an array via the Fischer-Yates algorithm.
     * See https://stackoverflow.com/questions/1519736/random-shuffling-of-an-array
     * </p>
     * 
     * @param array - the array to shuffle
     * @param random - the random-number generator
     */
    private static void shuffleArray(final long[] array, final Random rand)
    {
        int index;
        for (int i = array.length - 1; i > 0; i--)
        {
            index = rand.nextInt(i + 1);
            if (index != i)
            {
                array[index] ^= array[i];
                array[i] ^= array[index];
                array[index] ^= array[i];
            }
        }
    }

    /**
     * <p>
     * The {@code calMeanStd} method calculates the mean
     * and standard deviation from an array of samples.
     * </p>
     * 
     * @param mean - a double[] array for returning the mean by reference
     * @param random - a double[] array for returning the std. dev. by reference
     */
    private static void calcMeanStd(final double[] samples,
                                    final double[] mean,
                                    final double[] std)
    {
        double sum = 0, sum2 = 0;
        for (int i = 0; i < samples.length; ++i) {
            sum += samples[i];
            sum2 += samples[i] * samples[i];
        }
        double n = (double) samples.length;
        mean[0] = sum / n;
        std[0] = Math.sqrt((n * sum2) - (sum * sum)) / n;
    }

    /**
     * <p>
     * Define a simple data set then build a k-d tree.
     * </p>
     */
    public static void main(String[] args) {

        // Set the defaults then parse the input arguments.
        int iterations = 1;
        int numPoints = 262144;
        int extraPoints = 100;
        int numDimensions = 3;
        int numThreads = 1;
        long searchDivisor = 10;
        int maximumNumberOfNodesToPrint = 5;
        int numNearestNeighbors = 5;
        boolean neighbors = false;
        boolean region = false;
		
        for (int i = 0; i < args.length; i++) {
            if ( args[i].equals("-i") || args[i].equals("--iterations") ) {
                iterations = Integer.parseInt(args[++i]); 
                continue;
            }
             if ( args[i].equals("-n") || args[i].equals("--numPoints") ) {
                numPoints = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-x") || args[i].equals("--extraPoints") ) {
                extraPoints = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-d") || args[i].equals("--numDimensions") ) {
                numDimensions = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-t") || args[i].equals("--numThreads") ) {
                numThreads = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-s") || args[i].equals("--searchDivisor") ) {
                searchDivisor = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-p") || args[i].equals("--maximumNodesToPrint") ) {
                maximumNumberOfNodesToPrint = Integer.parseInt(args[++i]); 
                continue;
            }
            if ( args[i].equals("-m") || args[i].equals("--nearestNeighborCount") ) {
                numNearestNeighbors = Integer.parseInt(args[++i]);
                continue;
            }
            if (args[i].equals("-g") || args[i].equals("--neighbors")) {
            neighbors = !neighbors;
            continue;
            }
            if (args[i].equals("-j") || args[i].equals("--region")) {
            region = !region;
            continue;
            }
            if (args[i].equals("-h") || args[i].equals("--help")) {
                System.out.println("\nUsage:\n\n");
                System.out.println("java TestKdTree [-n N] [-x X] [-d D] [-t T] [-c C] " +
                                   "[-g] [-m M] [-j] [-s S] [-p P] [-i] [-h]\n\n" +
                                   "where the command-line options are interpreted as follows.\n\n");
                System.out.println("-n The number N of randomly generated points used to build the k-d tree\n\n");
                System.out.println("-x The number X of duplicate points added to to randomly generated points\n\n");
                System.out.println("-d The number of dimensions D (aka k) of the k-d tree\n\n");
                System.out.println("-t The number of threads T used to build and search the static k-d tree\n\n");
                System.out.println("-c The multi-thread cutoff below which only a single thread is used to build the tree\n\n");
                System.out.println("-g Find nearest neighbors to a query point\n\n");
                System.out.println("-m The maximum number M of nearest neighbors to be found\n\n");
                System.out.println("-j Perform a region search in a hypercube centered at a query point\n\n");
                System.out.println("-s The search divisor S that modifies the size of the hypercube for region search\n\n");
                System.out.println("-p The maximum number P of nodes to report when reporting region search results\n\n");
                System.out.println("-i The number I of iterations of k-d tree creation\n\n");
                System.out.println("-h List the command-line options\n\n");;
                System.exit(0);
            }
            throw new RuntimeException("illegal command-line argument: " + args[i]);
        }

        // Calculate a delta coordinate by dividing the positive range of long
        // by the number of points and truncating the quotient. Because the positive
        // range is less than half the full range of long, multiplying the
        // delta coordinate by the number of points ought to produce a product
        // that is less than half the full range of long and therefore avoid
        // possible overflow when comparing keys via the superKeyCompare function.
        // Calculate a padding coordinate to center the coordinates about zero.
        final long deltaCoordinate = Long.MAX_VALUE / numPoints;
        final long padCoordinate = (Long.MAX_VALUE - (numPoints * deltaCoordinate)) / 2;

        // Declare and initialize the oneCoordinate array. Equally space each coordinate
        // across the range of centered coordinates.
        final long beginCoordinate = Long.MIN_VALUE + padCoordinate;
        long thisCoordinate = beginCoordinate;
        long endCoordinate = 0;
        final long[] oneCoordinate = new long[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            oneCoordinate[i] = thisCoordinate;
            endCoordinate = thisCoordinate;
            thisCoordinate += deltaCoordinate;
        }

        // These two coordinates indicate whether integer overlow occurs for region search.
        long maxCoordinate = endCoordinate * 2;
        long minCoordinate = -beginCoordinate * 2;

        System.out.println("\ndeltaCoordinate = " + deltaCoordinate);
        System.out.println("padCoordinate = " + padCoordinate);
        System.out.println("beginCoordinate = " + beginCoordinate);
        System.out.println("endCoordinate = " + endCoordinate);
        System.out.println("minCoordinate = " + minCoordinate);
        System.out.println("maxCoordinate = " + maxCoordinate);

        // Calculate the number of child threads to be the number of threads minus 1, then
        // calculate the maximum tree depth at which to launch a child thread.  Truncate
        // this depth such that the total number of threads, including the master thread, is
        // an integer power of 2, hence simplifying the launching of child threads by restricting
        // them to only the < branch of the tree for some number of levels of the tree.
        int n = 0;
        if (numThreads > 0) {
            while (numThreads > 0) {
                n++;
                numThreads >>= 1;
            }
            numThreads = 1 << (n - 1);
        } else {
            numThreads = 0;
        }
        final int childThreads = numThreads - 1;
        int maximumSubmitDepth = -1;
        if (numThreads < 2) {
            maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
        } else if (numThreads == 2) {
            maximumSubmitDepth = 0;
        } else {
            maximumSubmitDepth = (int) Math.floor( Math.log( (double) childThreads ) / Math.log(2.) );
        }
        System.out.println("\nNumber of child threads = " + childThreads + "  maximum submit depth = " + maximumSubmitDepth + "\n");
		
        // Create a fixed thread pool ExecutorService.
        ExecutorService executor = null;
        if (childThreads > 0) {
            try {
                executor = Executors.newFixedThreadPool(childThreads);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("executor exception " + e.getMessage());
            }
        }
		
        // Allocate arrays to store the execution times and statistics.
        final double[] createTime = new double[iterations];
        final double[] initTime = new double[iterations];
        final double[] sortTime = new double[iterations];
        final double[] removeTime = new double[iterations];
        final double[] buildTime = new double[iterations];
        final double[] verifyTime = new double[iterations];
        final double[] searchTime = new double[iterations];
        final double[] containsTime = new double[iterations];
        final double[] neighborsSearchTime = new double[iterations];
        final double[] neighborsBruteTime = new double[iterations];
        final double[] regionSearchTime = new double[iterations];
        final double[] regionBruteTime = new double[iterations];

        // Initialize the random-number generator using Pi as a seed.
        Random rand = new Random();
        rand.setSeed(3141592653589793239L);
		
        // Declare and initialize the coordinates array. Each element of the
        // coordinates array is a pair wherein the key is an array of (x, y, z, w...)
        // coordinates and the value is the string representation of an integer.
        extraPoints = (extraPoints < numPoints) ? extraPoints : numPoints - 1;
        final Pair[] coordinates = new Pair[numPoints + extraPoints];
        for (int i = 0; i < coordinates.length; ++i) {
            final long[] key = new long[numDimensions];
            final String value = String.valueOf(i);
            coordinates[i] = new Pair(key, value);
        }

        // Create query arrays for searching the k-d tree via region and nearest-neighbors searches.
        long[] query = new long[numDimensions];
        long[] queryLower = new long[numDimensions];
        long[] queryUpper = new long[numDimensions];
        for (int i = 0; i < numDimensions; i++) {
            query[i] = i;
            queryLower[i] = query[i] + (beginCoordinate / searchDivisor);
            queryUpper[i] = query[i] + (endCoordinate / searchDivisor);
        }
        // Ensure that each query lower bound <= the corresponding query upper bound.
        for (int i = 0; i < queryLower.length; ++i) {
            if (queryLower[i] > queryUpper[i]) {
                long tmp = queryLower[i];
                queryLower[i] = queryUpper[i];
                queryUpper[i] = tmp;
            }
        }

        // These variables will be modified by building and testing the k-d tree.
        int treeHeight = 0, staticTreeHeight = 0, numberOfNodes = 0, staticNumberOfNodes;
        int numRegionNodes = 0, numNeighborsNodes = 0;

        // Build and test the k-d tree for the specified number of iterations.
        for (int k = 0; k < iterations; ++k)
        {
            // Shuffle the coordinates vector independently for each dimension.
            for (int j = 0; j < numDimensions; ++j) {
                shuffleArray(oneCoordinate, rand);
                for (int i = 0; i < numPoints; ++i) {
                    coordinates[i].getKey()[j] = oneCoordinate[i];
                }
            }

            // Reflect tuples across coordinates[numPoints - 1] to initialize the extra points.
            for (int i = 1; i <= extraPoints; ++i) {
                for (int j = 0; j < numDimensions; ++j) {
                    coordinates[numPoints - 1 + i].getKey()[j] = coordinates[numPoints - 1 - i].getKey()[j];
                }
            }

            // Create a static, balanced k-d tree.
            final long[] nN = new long[1];
            final double[] iT = new double[1]; // initialization time
            final double[] sT = new double[1]; // sort time
            final double[] rT = new double[1]; // remove time
            final double[] kT = new double[1]; // k-d tree-build time
            final double[] vT = new double[1]; // verify time

            final KdTree tree = KdTree.createKdTree(coordinates, executor, maximumSubmitDepth,
                                                    nN, iT, sT, rT, kT, vT);

            createTime[k] = iT[0] + sT[0] + rT[0] + kT[0] + vT[0];
            initTime[k] = iT[0];
            sortTime[k] = sT[0];
            buildTime[k] = kT[0];

            // Verify correct order of each node in the k-d tree and count the nodes.
            long vTime = System.currentTimeMillis();
            numberOfNodes = tree.verifyKdTree();
            vTime = System.currentTimeMillis() - vTime;
            verifyTime[k] += (double) vTime / Constants.MILLISECONDS_TO_SECONDS;

            if (numberOfNodes + extraPoints != coordinates.length) {
                throw new RuntimeException("number of coordinates = " + coordinates.length + "  !=  " +
                                           "number of nodes + extra points = " + (numberOfNodes + extraPoints));
            } else {
                treeHeight = tree.getTreeHeight();
            }

            // Check that the two versions of verifyKdTree report the same number of nodes.
            if (nN[0] != numberOfNodes) {
                throw new RuntimeException("number of nodes from createKdTree = " + nN[0] +
                                        "  != number of nodes from returned root = " + numberOfNodes);
            }
            
            // Search for each coordinate in the k-d tree. 
            long sTime = System.currentTimeMillis();
            for (int i = 0; i < coordinates.length; ++i) {
                if (!tree.contains(coordinates[i])) {
                    throw new RuntimeException("failed to find tuple " + i);
                }
            }
            sTime = System.currentTimeMillis() - sTime;
            searchTime[k] += (double) sTime / Constants.MILLISECONDS_TO_SECONDS;

            // Region search the k-d tree to find the nodes within a hypercube search region centered near the origin.
            if (region) {
                // Search the tree to get the list of KdNodes
                long rsTime = System.currentTimeMillis();
                List<KdNode> regionNodes = tree.searchKdTree(queryLower, queryUpper, executor, maximumSubmitDepth, 0, 0, true);
                rsTime = System.currentTimeMillis() - rsTime;
                regionSearchTime[k] += (double) rsTime / Constants.MILLISECONDS_TO_SECONDS;
                numRegionNodes = regionNodes.size();

                // Search the tree again to get the list of KdNodes.
                long bsTime = System.currentTimeMillis();
                List<KdNode> bruteNodes = tree.searchKdTree(queryLower, queryUpper, executor, maximumSubmitDepth, 0, 0, false);
                bsTime = System.currentTimeMillis() - bsTime;
                regionBruteTime[k] += (double) bsTime / Constants.MILLISECONDS_TO_SECONDS;

                // Compare the results of region search and brute-force search.
                if (regionNodes.size() != bruteNodes.size()) {
                    throw new RuntimeException("number of nodes found by region-search and brute-force do not match");
                } else {
                    for (int i = 0; i < regionNodes.size(); ++i) {
                        if (MergeSort.superKeyCompare(regionNodes.get(i).tuple, bruteNodes.get(i).tuple, 0) != 0L) {
                            throw new RuntimeException("region-search and brute-force values at " + i + " do not match");
                        }
                    }
                }
            }

            // Search the k-d tree to find the numNearestNeighbors nearest neighbors to the first point.
            if (neighbors)
            {
                // Search the tree to get the list of KdNodes.
                long nnTime = System.currentTimeMillis();
                NearestNeighborList nnList = new NearestNeighborList(query, numNearestNeighbors);
                tree.nearestNeighbor(nnList, 0);
                nnTime = System.currentTimeMillis() - nnTime;
                neighborsSearchTime[k] += (double) nnTime / Constants.MILLISECONDS_TO_SECONDS;
                numNeighborsNodes = nnList.curDepth;

                // Search the tree again to get the list of KdNodes.
                long bfTime = System.currentTimeMillis();
                NearestNeighborList bfList = new NearestNeighborList(query, numNearestNeighbors);
                tree.bruteNeighbor(bfList, 0);
                bfTime = System.currentTimeMillis() - bfTime;
                neighborsBruteTime[k] += (double) bfTime / Constants.MILLISECONDS_TO_SECONDS;

                // Compare the results of nearest-neighbor search and brute-force search.
                for (int i = 0; i < numNearestNeighbors; ++i) {
                    if (MergeSort.superKeyCompare(bfList.nodes[i].tuple, nnList.nodes[i].tuple, 0) != 0L) {
                        System.out.println("nearest-neighbor and brute-force values at " + i + " do not match");
                        System.out.println("nn dist = " + nnList.dists[i].toString() +
                                           "  bf dist = " + bfList.dists[i].toString() + "\n");
                    }
                }
            }

            System.out.println("finished iteration " + (k + 1));
        }

        // Report k-d tree statistics.
        System.out.println("\nstatic tree:\n");
        System.out.println("number of nodes = " + numberOfNodes + "  k-d tree height = " + treeHeight + "\n");

        // Here is a kludge for passing parameters to the calcMeanStd method by reference.
        final double[] mean = new double[1];
        final double[] std = new double[1];

        // Report execution times.
        calcMeanStd(createTime, mean, std);
        System.out.printf("create time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(verifyTime, mean, std);
        System.out.printf("verify time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(searchTime, mean, std);
        System.out.printf("search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(initTime, mean, std);
        System.out.printf("\ninit time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(sortTime, mean, std);
        System.out.printf("sort time = %.4f  std dev = %.4f\n", mean[0], std[0]);

        if (neighbors) {
            System.out.print("\nFound " + numNeighborsNodes + " nearest neighbors to ");
            KdNode.printTuple(query);
            System.out.println("\n");
            calcMeanStd(neighborsSearchTime, mean, std);
            System.out.printf("nearest neighbor search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
            calcMeanStd(neighborsBruteTime, mean, std);
            System.out.printf("nearest neighbor brute time  = %.4f  std dev = %.4f\n", mean[0], std[0]);
        }

        if (region) {
            System.out.print("\nFound " + numRegionNodes+ " tuples found by region search within " +
                              (queryUpper[0] - queryLower[0]) + " units of ");
            KdNode.printTuple(query);
            System.out.println(" in all dimensions.\n");
            calcMeanStd(regionSearchTime, mean, std);
            System.out.printf("region search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
            calcMeanStd(regionBruteTime, mean, std);
            System.out.printf("brute search time  = %.4f  std dev = %.4f\n", mean[0], std[0]);
        }

        System.out.println();
    
        // Shut down the ExecutorService.
        if (childThreads > 0) {
            executor.shutdown();
        }
    }
} // class TestKdTree
