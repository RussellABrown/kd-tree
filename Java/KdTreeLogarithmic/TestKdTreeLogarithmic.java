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

/*
 * Test program for KdTreeLogarithmic.java, KdTreeDynamic.java, KdTree.java, KdTreeNlogn.java,
 * KdTreeKnlogn.java, KdNode.java, MergeSort.java, NearestNeighborHeap.java, Pair.java, Paire.java,
 * AvlTree.java, AvlNode.java and Constants.java
 *
 * Configuration is controlled via the following constants in Constants.java
 * 
 * NLOGN - If true, specifies the O[n log(n)] static k-d tree algorithm instead of
 *         the O[kn log(n)] static k-d tree algorithm.
 * 
 * INSERTION_SORT_CUTOFF - An integer that specifies the number of nodes at which insertion
 *                         sort replaces merge sort in the MergeSort.mergeSort* methods.
 * 
 * MEDIAN_OF_MEDIANS_CUTOFF - An integer that specifies number of nodes at which insertion
 *                            sort replaces the median-of-medians algorithm in the
 *                            KdTreeNlogn.partiiton method.
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
 * MULTI_THREAD_CUTOFF - An integer that specifies the number of nodes above which
 *                       multiple threads are used to construct a static k-d tree.
 * 
 * AVL_BALANCE - If true, the KdTreeDynamic.isBalanced method checks for AVL balance;
 *               otherwise, this method checks for red-black balance. It appears that
 *               red-black balance confers better performance than AVL balance.
 * 
 * HEIGHT_DIFF - An integer that specifies, for red-black balance, the maximum allowed
 *               height difference between the < and > subtrees of a node when one subtree
 *               is empty; or for AVL balance, the maximum allowed height difference
 *               between the < and > subtrees of a node. This balance is used to maintain
 *               the balance of a dynamic k-d tree, not a static k-d tree.
 * 
 * ENABLE_DEBUG - If true, perform additional validity checks.
 * 
 * ENABLE_1TO3 - If true, rebalancing of a dynamic k-d tree is performed using simple
 *               operations, and recursive k-d deletion is curtailed, when a subtree
 *               contains <= 3 nodes.
 * 
 * ENABLE_PREFERRED_AVL_NODE - If true, inspect the balance of a deleted 2-child AVL node
 *                             to select a preferred replacement node.
 * 
 * ENABLE_PREFERRED_KD_NODE - If true, compare the heights of a deleted 2-child k-d node's
 *                            child subtrees to select a preferred replacement node.
 * 
 * ENABLE_LINKED_LIST - If true, use the LinkedList version of region search
 *                      instead of the ArrayList version.
 * 
 * ENABLE_INSERTION_REBALANCE - If true, a dynamic k-d tree is rebalanced after insertion
 *                              of a node if necessary.
 * 
 * ENABLE_DELETION_REBALANCE - If true, a dynamic k-d tree is rebalanced after deletion
 *                             of a node if necessary.
 * 
 * ENABLE_HISTOGRAMS - If true, histograms are obtained for rebalancing operations.
 * 
 * ENABLE_SPARSE_INSERTION - If true, the KdTreeDynamic.insert method inserts a node
 *                           into the smallest non-full dynamic k-d tree of the set
 *                           of k-d trees managed by a logarithmic k-d tree.
 * 
 * ENABLE_TUPLE_COPY - If true, specifies that the tuple array is copied in the
 *                     KdTreeDynamic.erase method and the KdNode and AvlNode constructors.
 * 
 * ENABLE_LIST_PREPEND - If true, a new k-d node is prepended to the doubly linked list
 *                       by the KdTree.add mthod; otherwise, the new k-d node is appended.
 * 
 * ENABLE_LIST_PREPEND_ALL - If true, a doubly linked list of k-d nodes is prepended to
 *                           another doubly linked list by KdTree.addList method; otherwise,
 *                           the list is appended.
 * 
 * 
 *
 * 
 * Usage:
 *
 * java TestKdTreeLogarithmic [-n N] [-x X] [-d D] [-t T] [-b] [-g] [-m M] [-j] [-s S] \
 *                            [-p P] [-v] [-f] [-c] [-w] [-y] [-z Z] [-r] [-i] [-h]
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
 * -v Verify a correct k-d tree after each insertion or erasure of each pair (default off)
 * 
 * -f Search for the next pair after erasure of each pair (default off)
 * 
 * -c Search for all remaining pairs after erasure of each pair (default off)
 * 
 * -w Create a worst-case set of points by walking a k-d tree in order (default off)
 * 
 * -y Create a set of points that lie on the diagonal of the hypercube (default off)
 * 
 * -z The number Z of segments into which to break the points for erasure/insertion (default 1)
 * 
 * -r Reverse the order of points for erasure relative to insertion (default off)
 *
 * -i The number I of iterations of k-d tree insertion, search, and deletion (default 1)
 *
 * -h Help, i.e., print the above explanation of command-line options.
 */

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
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
public class TestKdTreeLogarithmic {

    /**
     * <p>
     * The {@code shuffleArray} method shuffles an array via the Fischer-Yates algorithm.
     * See https://stackoverflow.com/questions/1519736/random-shuffling-of-an-array
     * </p>
     * 
     * @param array - the array to shuffle
     * @param random - the random-number generator
     */
    private static void shuffleArray(final long[] array,
                                     final Random rand)
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
    public static void main(String[] args)
    {
        // Set the defaults then parse the input arguments.
        int iterations = 1;
        int numPoints = 262144;
        int extraPoints = 100;
        int numDimensions = 3;
        int numThreads = 1;
        long searchDivisor = 10;
        int maximumNumberOfNodesToPrint = 5;
        int numNearestNeighbors = 5;
        int fraction = 1;
        boolean neighbors = false;
        boolean region = false;
        boolean balanced = false;
        boolean verify = false;
        boolean find = false;
        boolean check = false;
        boolean worst = false;
        boolean diagonal = false;
        boolean reverse = false;
		
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
            if (args[i].equals("-b") || args[i].equals("--balanced")) {
                balanced = !balanced;
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
            if (args[i].equals("-v") || args[i].equals("--verify")) {
                verify = !verify;
                continue;
            }
            if (args[i].equals("-f") || args[i].equals("--find")) {
                find = !find;
                continue;
            }
            if (args[i].equals("-c") || args[i].equals("--check")) {
                check = !check;
                continue;
            }
            if (args[i].equals("-w") || args[i].equals("--worst")) {
                worst = !worst;
                continue;
            }
            if (args[i].equals("-y") || args[i].equals("--diagonal")) {
                diagonal = !diagonal;
                continue;
            }
            if ( args[i].equals("-z") || args[i].equals("--fraction") ) {
                fraction= Integer.parseInt(args[++i]);
                continue;
            }
            if (args[i].equals("-r") || args[i].equals("--reverse")) {
                reverse = !reverse;
                continue;
            }
            if (args[i].equals("-h") || args[i].equals("--help")) {
                System.out.println("\nUsage:\n");
                System.out.println("java TestKdTreeLogarithmic [-n N] [-x X] [-d D] [-t T] [-b] [-g] [-m M] [-j] " +
                                   "[-s S] [-p P] [-v] [-f] [-c] [-w] [-y] [-z Z] [-r] [-i] [-h]\n\n" +
                                   "where the command-line options are interpreted as follows.\n");
                System.out.println("-n The number N of randomly generated points used to build the k-d tree\n");
                System.out.println("-x The number X of duplicate points added to to randomly generated points\n");
                System.out.println("-d The number of dimensions D (aka k) of the k-d tree\n");
                System.out.println("-t The number of threads T used to build and search the static k-d tree\n");
                System.out.println("-b Build a balanced, static k-d tree for comparison to the dynamic k-d tree\n");
                System.out.println("-g Find nearest neighbors to a query point\n");
                System.out.println("-m The maximum number M of nearest neighbors to be found\n");
                System.out.println("-j Perform a region search in a hypercube centered at a query point\n");
                System.out.println("-s The search divisor S that modifies the size of the hypercube for region search\n");
                System.out.println("-p The maximum number P of nodes to report when reporting region search results\n");
                System.out.println("-v Verify the k-d tree ordering and balance after insertion or erasure of each pair\n");
                System.out.println("-f Check for the next point after deleting each point (a cheap alternative to -v)\n");
                System.out.println("-c Check for all remaining points after deleting each point (an expensive alternative to -v)\n");
                System.out.println("-w Create a worst-case set of point by walking a k-d tree in order\n");
                System.out.println("-y Create a set of points that lie on the diagonal of the hypercube\n");
                System.out.println("-z The number Z of segments into which to break the points for erasure/insertion\n");
                System.out.println("-r Reverse the order of points for erasure relative to insertion\n");
                System.out.println("-i The number I of iterations of k-d tree creation\n");
                System.out.println("-h List the command-line options\n");
                System.exit(0);
            }
            throw new RuntimeException("\n\nillegal command-line argument: " + args[i] + "\n");
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
        final double[] insertTime = new double[iterations];
        final double[] verifyTime = new double[iterations];
        final double[] searchTime = new double[iterations];
        final double[] eraseTime = new double[iterations];
        final double[] finsertTime = new double[iterations];
        final double[] feraseTime = new double[iterations];
        final double[] containsTime = new double[iterations];
        final double[] neighborsSearchTime = new double[iterations];
        final double[] neighborsBruteTime = new double[iterations];
        final double[] regionSearchTime = new double[iterations];
        final double[] regionBruteTime = new double[iterations];

        // Initialize the random-number generator using Pi as a seed. An alternative to Random might be the Mersenne twister:
        // https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/MersenneTwister.html
        Random rand = new Random();
        rand.setSeed(3141592653589793239L);
		
        // Declare and initialize the coordinates array. Each element of the
        // coordinates array is a pair wherein the key is an array of (x, y, z, etc.)
        // coordinates and the value is the string representation of an integer.
        //
        // Do not append extra pairs if worst-case (sorted) order is specified,
        // because worst-case order is achieved by building a k-d tree from
        // randomly ordered pairs and then walking the tree in sorted order
        // in order to cull sorted-order pairs. Because extra pairs have
        // duplicate tuples their values are stored in a key-to-multiple-values
        // map. However, culling sorted-order pairs takes only the first
        // of multiple values to form the pair. Hence, not all values will
        // be represented in the sorted-order pairs. This problem is avoided
        // by appending no extra pairs when worst-case order is specified.
        extraPoints = (extraPoints < numPoints) ? extraPoints : numPoints - 1;
        if (worst) {
            extraPoints = 0;
        }
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

        // Create an instance of KdTreeLogarithmic.
        final KdTreeLogarithmic tree = new KdTreeLogarithmic(numDimensions, executor, maximumSubmitDepth);

        // These variables will be modified by building and testing the k-d tree.
        long treeSize = 0L, numberOfNodes = 0L, staticNumberOfNodes = 0L;
        int treeHeight = 0, staticTreeHeight = 0, numRegionNodes = 0, numNeighborsNodes = 0;

        // Build and test the k-d tree for the specified number of iterations.
        for (int k = 0; k < iterations; ++k)
        {
            // Shuffle the coordinates vector independently for each dimension,
            // unless a diagonal set of coordinates is specified.
            if (diagonal) {
                shuffleArray(oneCoordinate, rand);
                for (int j = 0; j < numDimensions; ++j) {
                    for (int i = 0; i < numPoints; ++i) {
                        coordinates[i].getKey()[j] = oneCoordinate[i];
                    }
                }
            } else {
                for (int j = 0; j < numDimensions; ++j) {
                    shuffleArray(oneCoordinate, rand);
                    for (int i = 0; i < numPoints; ++i) {
                        coordinates[i].getKey()[j] = oneCoordinate[i];
                    }
                }
            }

            // Reflect pairs across coordinates[numPoints - 1] to initialize the extra pairs.
            for (int i = 1; i <= extraPoints; ++i) {
                for (int j = 0; j < numDimensions; ++j) {
                    coordinates[numPoints - 1 + i].getKey()[j] = coordinates[numPoints - 1 - i].getKey()[j];
                }
            }

            // Wrap a static k-d tree in a dynamic k-d tree. Walk the tree
            // in sorted order and over-write the first numPoints pairs
            // in the coordinates array to create worst-case coordinates.
            // Return the number ofnodes and execution times by reference
            // via single-element arrays.
            if (worst) {
                final long[] nN = new long[1];
                final double[] iT = new double[1]; // initialization time
                final double[] sT = new double[1]; // sort time
                final double[] rT = new double[1]; // remove time
                final double[] kT = new double[1]; // k-d tree-build time
                final double[] vT = new double[1]; // verify time

                final KdTree arbre = new KdTreeDynamic(numDimensions, executor,
                                                       maximumSubmitDepth,
                                                       KdTree.createKdTree(coordinates,
                                                                           executor,
                                                                           maximumSubmitDepth,
                                                                           nN,
                                                                           iT,
                                                                           sT,
                                                                           rT,
                                                                           kT,
                                                                           vT));

                // Check that the static tree contains the correct number of nodes.
                long worstNumberOfNodes = arbre.verifyKdTree();
                if (nN[0] != worstNumberOfNodes) {
                    throw new RuntimeException("\n\nfor worst tree, number of nodes from createKdTree = " + nN[0] +
                                               "  != number of nodes counted by verifyKdTree = " +
                                               worstNumberOfNodes + "\n");
                }
                if (numPoints != worstNumberOfNodes) {
                    throw new RuntimeException("\n\nfor worst tree, number of non-extra pairs = " + nN[0] +
                                               "  != number of nodes counted by verifyKdTree = " +
                                               worstNumberOfNodes + "\n");
                }


                // Walk the tree in sorted order and over-write the first
                // numPoints pairs in the coordinates array. Pass the 'index'
                // counter by referenced via a single-element array.
                int[] index = {0};
                int worstPoints = tree.getSortedPairs(arbre.root, coordinates, index);

                if (worstPoints != index[0]) {
                    throw new RuntimeException("\n\nnumber of sorted points = " + worstPoints + "  !=  " +
                                               "number of indexed points = " + index[0] + "\n");
                }

                if (worstPoints != numPoints) {
                    throw new RuntimeException("\n\nnumber of non-extra coordinates = " + numPoints + "  !=  " +
                                               "number of sorted points = " + worstPoints + "\n");    
                }
            }

           // Insert each coordinate into the dynamic k-d tree.
            long iTime = System.currentTimeMillis();
            for (int i = 0; i < coordinates.length; ++i) {
                if (tree.insert(coordinates[i])) {
                    if (verify) {
                        tree.verifyKdTree();
                    }
                } else {
                    throw new RuntimeException("\n\nfailed to insert pair " + i + "\n");
                }
            }

            iTime = System.currentTimeMillis() - iTime;
            insertTime[k] += (double) iTime / Constants.MILLISECONDS_TO_SECONDS;

            // Verify correct order of each node in the k-d tree and count the nodes.
            long vTime = System.currentTimeMillis();
            numberOfNodes = tree.verifyKdTree();
            vTime = System.currentTimeMillis() - vTime;
            verifyTime[k] += (double) vTime / Constants.MILLISECONDS_TO_SECONDS;

            if (numberOfNodes + extraPoints != coordinates.length) {
                throw new RuntimeException("\n\nnumber of coordinates = " + coordinates.length + "  !=  " +
                                           "number of nodes + extra points = " + (numberOfNodes + extraPoints) + "\n");
            }

            // Record the maximum tree height and the tree size.
            treeHeight = tree.getTreeHeight();
            treeSize = KdTreeLogarithmic.getSize(tree);

            // Search for each coordinate in the k-d tree. 
            long sTime = System.currentTimeMillis();
            for (int i = 0; i < coordinates.length; ++i) {
                if (!tree.contains(coordinates[i])) {
                    throw new RuntimeException("\n\nfailed to find pair " + i + "\n");
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
                    throw new RuntimeException("\n\nnumber of nodes found by region-search and brute-force do not match\n");
                } else {
                    for (int i = 0; i < regionNodes.size(); ++i) {
                        if (MergeSort.superKeyCompare(regionNodes.get(i).tuple, bruteNodes.get(i).tuple, 0) != 0L) {
                            throw new RuntimeException("\n\nregion-search and brute-force values at " + i + " do not match\n");
                        }
                    }
                }
            }

            // Search the k-d tree to find the numNearestNeighbors nearest neighbors to the first point.
            if (neighbors)
            {
                // Search the tree to get the list of KdNodes.
                long nnTime = System.currentTimeMillis();
                List<Paire> nnList = tree.findNearestNeighbors(query, numNearestNeighbors);
                nnTime = System.currentTimeMillis() - nnTime;
                neighborsSearchTime[k] += (double) nnTime / Constants.MILLISECONDS_TO_SECONDS;
                numNeighborsNodes = nnList.size();

                // Search the tree again to get the list of KdNodes.
                long bfTime = System.currentTimeMillis();
                List<Paire> bfList = tree.findBruteNeighbors(query, numNearestNeighbors);
                bfTime = System.currentTimeMillis() - bfTime;
                neighborsBruteTime[k] += (double) bfTime / Constants.MILLISECONDS_TO_SECONDS;
                int numBruteNodes = bfList.size();

                // Compare the results of nearest-neighbor search and brute-force search.
                if (numNeighborsNodes != numBruteNodes) {
                    System.out.println("nearest-neighbor size = " + numNeighborsNodes +
                                       "  !=  brute-force size = " + numBruteNodes);
                }
                for (int i = 0; i < numNeighborsNodes; ++i) {
                    if (MergeSort.superKeyCompare(bfList.get(i).getValue().tuple, nnList.get(i).getValue().tuple, 0) != 0L) {
                        System.out.println("nearest-neighbor and brute-force values at " + i + " do not match");
                        System.out.println("nn dist = " + nnList.get(i).getKey() +
                                           "  bf dist = " + bfList.get(i).getKey()+ "\n");
                    }
                }
            }

            // Erase and re-insert the coordinates in pieces.
            if (fraction > 1) {
                int size = coordinates.length / fraction;

                // Each piece must contain at least 1 coordinate.
                if (size >= 1) {

                    // Randomly shuffle the coordinates.
                    List<Pair> coordinateList = Arrays.asList(coordinates);
                    Collections.shuffle(coordinateList);
                    Pair[] shuffledCoordinates = coordinateList.toArray(new Pair[0]);

                    // Iterate over the pieces.
                    for (int i = 0, start = 0; i < fraction; ++i, start += size) {
                        // Prevent the end of the piece from exceeding the number of shuffled coordinates.
                        int end = (start + size > shuffledCoordinates.length) ? shuffledCoordinates.length : start + size;
                        if (reverse) {
                            long feraTime = System.currentTimeMillis();
                            for (int j = end-1; j >= start; --j) {
                                if (tree.erase(shuffledCoordinates[j])) {
                                    // Verify correctness of the k-d tree after each erasure.
                                    if (verify) {
                                        tree.verifyKdTree();
                                    }
                                    // A search for a (key, value) pair after erasing it should fail.
                                    if (find && tree.contains(shuffledCoordinates[j])) {
                                        throw new RuntimeException("\n\nfound fractional pair " + j + " after erasing it\n");
                                    }
                                    // A search for the prior (key, value) pair after erasing a (key, value) pair should succeed.
                                    if (find && j > 0 && !tree.contains(shuffledCoordinates[j-1])) {
                                        throw new RuntimeException("\n\nfailed to find prior fractional pair " + (j-1) + " after erasing pair " + j + "\n");
                                    }
                                    // A search for all prior (key value) pairs after erasing a (key, value) pair should succeed.
                                    if (check && j > 0) {
                                        for (int l = j-1; l >= 0; --l) {
                                            if (!tree.contains(coordinates[l])) {
                                                throw new RuntimeException("\n\nfailed to find prior fractional pair " + l + " after erasing pair " + j + "\n");
                                            }
                                        }
                                    }
                                } else {
                                    throw new RuntimeException("\n\nfailed to erase fractional pair " + j + " from dynamic k-d tree\n");
                                }
                            }
                            feraTime = System.currentTimeMillis() - feraTime;
                            feraseTime[k] += (double) feraTime / Constants.MILLISECONDS_TO_SECONDS;
                            long finsTime = System.currentTimeMillis();
                            for (int j = start; j < end; ++j) {
                                if (tree.insert(shuffledCoordinates[j])) {
                                    if (verify) {
                                        tree.verifyKdTree();
                                    }
                                } else {
                                    throw new RuntimeException("\n\nfailed to insert fractional pair " + j + "\n");
                                }
                            }
                            finsTime = System.currentTimeMillis() - finsTime;
                            finsertTime[k] += (double) finsTime / Constants.MILLISECONDS_TO_SECONDS;
                        } else {
                            long feraTime = System.currentTimeMillis();
                            for (int j = start; j < end; ++j) {
                                if (tree.erase(shuffledCoordinates[j])) {
                                    // Verify correctness of the k-d tree after each erasure.
                                    if (verify) {
                                        tree.verifyKdTree();
                                    }
                                    // A search for a (key, value) pair after erasing it should fail.
                                    if (find && tree.contains(shuffledCoordinates[j])) {
                                        throw new RuntimeException("\n\nfound fractional pair " + j + " after erasing it\n");
                                    }
                                    // A search for the next (key, value) pair after erasing a (key, value) pair should succeed.
                                    if (find && j < coordinates.length-1 && !tree.contains(coordinates[j+1])) {
                                        throw new RuntimeException("\n\nfailed to find next fractional pair " + (j+1) + " after erasing pair " + j + "\n");
                                    }
                                    // A search for all following (key, value) pairs after erasing a (key, value) pair should succeed.
                                    if (check && j < coordinates.length-1) {
                                        for (int l = j+1; l < coordinates.length; ++l) {
                                            if (!tree.contains(coordinates[l])) {
                                                throw new RuntimeException("\n\nfailed to find following fractional pair " + l + " after erasing pair " + j + "\n");
                                            }
                                        }
                                    }
                                } else {
                                    throw new RuntimeException("\n\nfailed to erase fractional pair " + j + " from dynamic k-d tree\n");
                                }
                            }
                            feraTime = System.currentTimeMillis() - feraTime;
                            feraseTime[k] += (double) feraTime / Constants.MILLISECONDS_TO_SECONDS;
                            long finsTime = System.currentTimeMillis();
                            for (int j = start; j < end; ++j) {
                                if (tree.insert(shuffledCoordinates[j])) {
                                    if (verify) {
                                        tree.verifyKdTree();
                                    }
                                } else {
                                    throw new RuntimeException("\n\nfailed to insert fractional pair " + j + "\n");
                                }
                            }
                            finsTime = System.currentTimeMillis() - finsTime;
                            finsertTime[k] += (double) finsTime / Constants.MILLISECONDS_TO_SECONDS;
                        }
                    }
                }
            }

            // Erase each coordinate from the dynamic k-d tree, and reverse
            // the order of the coordinates if reverse is true.
            long eTime = System.currentTimeMillis();
            if (reverse) {
                for (int i = coordinates.length - 1; i >= 0; --i) {
                    if (tree.erase(coordinates[i])) {
                        // Verify correctness of the k-d tree after each erasure.
                        if (verify) {
                            tree.verifyKdTree();
                        }
                        // A search for a (key, value) pair after erasing it should fail.
                        if (find && tree.contains(coordinates[i])) {
                            throw new RuntimeException("\n\nfound pair " + i + " after erasing it\n");
                        }
                        // A search for the prior (key, value) pair after erasing a (key, value) pair should succeed.
                        if (find && i > 0 && !tree.contains(coordinates[i-1])) {
                            throw new RuntimeException("\n\nfailed to find prior pair " + (i-1) + " after erasing pair " + i + "\n");
                        }
                        // A search for all prior (key value) pairs after erasing a (key, value) pair should succeed.
                        if (check && i > 0) {
                            for (int j = i-1; j >= 0; --j) {
                                if (!tree.contains(coordinates[j])) {
                                    throw new RuntimeException("\n\nfailed to find prior pair " + j + " after erasing pair " + i + "\n");
                                }
                            }
                        }
                    } else {
                        throw new RuntimeException("\n\nfailed to erase pair " + i + " from dynamic k-d tree\n");
                    }
                }
            } else {
                for (int i = 0; i < coordinates.length; ++i) {
                    if (tree.erase(coordinates[i])) {
                        // Verify correctness of the k-d tree after each erasure.
                        if (verify) {
                            tree.verifyKdTree();
                        }
                        // A search for a (key, value) pair after erasing it should fail.
                        if (find && tree.contains(coordinates[i])) {
                            throw new RuntimeException("\n\nfound pair " + i + " after erasing it\n");
                        }
                        // A search for the next (key, value) pair after erasing a (key, value) pair should succeed.
                        if (find && i < coordinates.length-1 && !tree.contains(coordinates[i+1])) {
                            throw new RuntimeException("\n\nfailed to find next pair " + (i+1) + " after erasing pair " + i + "\n");
                        }
                        // A search for all following (key, value) pairs after erasing a (key, value) pair should succeed.
                        if (check && i < coordinates.length-1) {
                            for (int j = i+1; j < coordinates.length; ++j) {
                                if (!tree.contains(coordinates[j])) {
                                    throw new RuntimeException("\n\nfailed to find following pair " + j + " after erasing pair " + i + "\n");
                                }
                            }
                        }
                    } else {
                        throw new RuntimeException("\n\nfailed to erase pair " + i + " from dynamic k-d tree\n");
                    }
                }
            }
            eTime = System.currentTimeMillis() - eTime;
            eraseTime[k] += (double) eTime / Constants.MILLISECONDS_TO_SECONDS;

            if ( !tree.isEmpty() ) {
                throw new RuntimeException("\n\ntree is not empty following erasure of all pairs\n");
            }

            // Wrap a static k-d tree in a dynamic k-d tree. Return the number of
            // nodes and execution times by reference via single-element arrays.
            if (balanced)
            {
                final long[] nN = new long[1];
                final double[] iT = new double[1]; // initialization time
                final double[] sT = new double[1]; // sort time
                final double[] rT = new double[1]; // remove time
                final double[] kT = new double[1]; // k-d tree-build time
                final double[] vT = new double[1]; // verify time

                final KdTree arbre = new KdTreeDynamic(numDimensions, executor,
                                                       maximumSubmitDepth,
                                                       KdTree.createKdTree(coordinates,
                                                                           executor,
                                                                           maximumSubmitDepth,
                                                                           nN,
                                                                           iT,
                                                                           sT,
                                                                           rT,
                                                                           kT,
                                                                           vT));

                staticTreeHeight = arbre.getTreeHeight();
                createTime[k] = iT[0] + sT[0] + rT[0] + kT[0];

                // Check that the static tree contains the correct number of nodes.
                staticNumberOfNodes = arbre.verifyKdTree();
                if (nN[0] != staticNumberOfNodes) {
                    throw new RuntimeException("\n\nfor static tree, number of nodes from createKdTree = " + nN[0] +
                                               "  != number of nodes counted by verifyKdTree = " +
                                               staticNumberOfNodes + "\n");
                }
                if (numPoints != staticNumberOfNodes) {
                    throw new RuntimeException("\n\nfor static tree, number of non-extra pairs = " + nN[0] +
                                               "  != number of nodes counted by verifyKdTree = " +
                                               staticNumberOfNodes + "\n");
                }

                // Search the dynamic k-d tree for each of the coordinates. Creation of
                // this dynamic k-d tree provides an example of how to provide a static
                // k-d tree to the dynamic k-d tree constructor, and thereafter treat
                // the static k-d tree as a dynamic k-d tree for various purposes, such
                // as insertion, deletion, and search of individual coordinates.
                long cTime = System.currentTimeMillis();
                for (int i = 0; i < coordinates.length; ++i) {
                        if (!arbre.contains(coordinates[i])) {
                            throw new RuntimeException("\n\nfailed to find pair " + i + " in static tree\n");
                        }
                    }
                cTime = System.currentTimeMillis() - cTime;
                containsTime[k] += (double) cTime / Constants.MILLISECONDS_TO_SECONDS;
            }

            System.out.println("finished iteration " + (k + 1));
        }

        // Report k-d tree statistics.
        System.out.println("\nlogarithmic tree:\n");
        System.out.println("total number of nodes = " + treeSize + "  maximum k-d tree height = " +
                           treeHeight + "\n");

        // Here is a kludge for passing parameters to the calcMeanStd method by reference.
        final double[] mean = new double[1];
        final double[] std = new double[1];

        // Report execution times.
        calcMeanStd(insertTime, mean, std);
        System.out.printf("insert time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(verifyTime, mean, std);
        System.out.printf("verify time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(searchTime, mean, std);
        System.out.printf("search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        calcMeanStd(eraseTime, mean, std);
        System.out.printf("delete time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        if (fraction > 1) {
            calcMeanStd(finsertTime, mean, std);
            System.out.printf("finser time = %.4f  std dev = %.4f\n", mean[0], std[0]);
            calcMeanStd(feraseTime, mean, std);
            System.out.printf("fdelet time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        }

        if (balanced) {
            System.out.println("\nstatic tree:\n");
            System.out.println("number of nodes = " + numberOfNodes + "  k-d tree height = " +
                               staticTreeHeight + "\n");
            calcMeanStd(createTime, mean, std);
            System.out.printf("create time = %.4f  std dev = %.4f\n", mean[0], std[0]);
            calcMeanStd(containsTime, mean, std);
            System.out.printf("search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
        }

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
            System.out.print("\nFound " + numRegionNodes + " pairs by region search within " +
                              (queryUpper[0] - queryLower[0]) + " units of ");
            KdNode.printTuple(query);
            System.out.println(" in all dimensions.\n");
            calcMeanStd(regionSearchTime, mean, std);
            System.out.printf("region search time = %.4f  std dev = %.4f\n", mean[0], std[0]);
            calcMeanStd(regionBruteTime, mean, std);
            System.out.printf("brute search time  = %.4f  std dev = %.4f\n", mean[0], std[0]);
        }

        if (Constants.ENABLE_HISTOGRAMS) {
            tree.sumHistograms();
            System.out.println("\nHistograms of built tree sizes follow.");
            System.out.println("\npow\tinsert log\tdelete log\tinsert dyn\tdelete dyn\tswap log\n");
            for (int i = 0; i < Constants.MAX_POWER_OF_2; ++i) {
                System.out.println(i + "\t" + tree.insertionHistogramLog[i] +
                                   "\t\t" + tree.deletionHistogramLog[i] +
                                   "\t\t" + tree.insertHistogramDyn[i] +
                                   "\t\t" + tree.deleteHistogramDyn[i] +
                                   "\t\t" + tree.swapHistogramLog[i]);
            }
            System.out.println("\nHistograms of n*log(n) operations follow.");
            System.out.println("\npow\tinsert log\tdelete log\tinsert dyn\tdelete dyn\n");
            System.out.println("0\t-Inf\t\t-Inf\t\t-Inf\t\t-Inf");
            int powerOf2 = 2;
            double insertSumLog = 0, deleteSumLog = 0, insertSumDyn = 0, deleteSumDyn = 0;
            for (int i = 1; i < Constants.MAX_POWER_OF_2; ++i) {
                double insertValueLog = (double) (tree.insertionHistogramLog[i]) * (double) (i * powerOf2);
                double deleteValueLog = (double) (tree.deletionHistogramLog[i]) * (double) (i * powerOf2);
                double insertValueDyn = (double) (tree.insertHistogramDyn[i]) * (double) (i * powerOf2);
                double deleteValueDyn = (double) (tree.deleteHistogramDyn[i]) * (double) (i * powerOf2);
                insertSumLog += insertValueLog;
                deleteSumLog += deleteValueLog;
                insertSumDyn += insertValueDyn;
                deleteSumDyn += deleteValueDyn;
                System.out.printf("%d\t%.2e\t%.2e\t%.2e\t%.2e\n", i,
                                  insertValueLog, deleteValueLog, insertValueDyn, deleteValueDyn);
            }
            System.out.println();
            System.out.println("total\tinsert log\tdelete log\tinsert dyn\tdelete dyn\n");
            System.out.printf("\t%.2e\t%.2e\t%.2e\t%.2e\n",
                              insertSumLog, deleteSumLog, insertSumDyn, deleteSumDyn);
        }
    
        System.out.println();

        // Shut down the ExecutorService.
        if (childThreads > 0) {
            executor.shutdown();
        }
    }
} // class TestKdTreeDynamic
