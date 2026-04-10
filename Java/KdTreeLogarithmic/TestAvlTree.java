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
 * Test program for AvlTree.java, Pair.java, Paire.java and Constants.java
 *
 * Configuration is controlled via the following constant in Constants.java
 * 
 * ENABLE_TUPLE_COPY - If true, specifies that the tuple array is copied in
 *                     the AvlNode constructor.
 *
 * 
 * Usage:
 *
 * "java TestAvlTree [-n N] [-x X] [-d D] [-v] [-f] [-r] [-i] [-h]
 *
 * where the command-line options are interpreted as follows.
 * 
 * -n The number N of randomly generated points used to build the AVL tree (default 262144)
 *
 * -x The number X of duplicate points added to the set of randomly generated points (default 100)
 *
 * -d The number of dimensions D (aka k) of the k-d tree (default 3)
 *
 * -v Verify a correct k-d tree after each insertion or erasure (default off)
 * 
 * -f Search for the tuple and the next tuple after erasure of each tuple (default off)
 * 
 * -r Reverse the order of coordinates for erasure relative to insertion (default off)
 *
 * -i The number I of iterations of k-d tree insertion, search, and deletion (default 1)
 *
 * -h Help, i.e., print the above explanation of command-line options.
 */

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
public class TestAvlTree {

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
        boolean verify = false;
        boolean find = false;
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
            if (args[i].equals("-v") || args[i].equals("--verify")) {
            verify = !verify;
            continue;
            }
            if (args[i].equals("-f") || args[i].equals("--find")) {
            find = !find;
            continue;
            }
            if (args[i].equals("-r") || args[i].equals("--reverse")) {
            find = !find;
            continue;
            }
            if (args[i].equals("-h") || args[i].equals("--help")) {
                System.out.println("\nUsage:\n");
                System.out.println("java TestKdTreeDynamic [-n N] [-x X] [-d D] [-v] [-f] [-r] [-i] [-h]\n\n" +
                                   "where the command-line options are interpreted as follows.\n");
                System.out.println("-n The number N of randomly generated points used to build the k-d tree\n");
                System.out.println("-x The number X of duplicate points added to to randomly generated points\n");
                System.out.println("-d The number of dimensions D (aka k) of the k-d tree\n");
                System.out.println("-v Verify the k-d tree ordering and balance after insertion or erasure of each point\n");
                System.out.println("-f Check for the next point after deleting each point (a cheap alternative to -v)\n");
                System.out.println("-r Reverse the order of coordinates for erasure relative to insertion\n");
                System.out.println("-i The number I of iterations of k-d tree creation\n");
                System.out.println("-h List the command-line options\n");
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
        System.out.println();

        // Allocate arrays to store the execution times and statistics.
        final double[] createTime = new double[iterations];
        final double[] insertTime = new double[iterations];
        final double[] verifyTime = new double[iterations];
        final double[] searchTime = new double[iterations];
        final double[] eraseTime = new double[iterations];

        // Initialize the random-number generator using Pi as a seed. An alternative to Random might be the Mersenne twister:
        // https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/MersenneTwister.html
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

        // Create an instance of AvlTree.
        final AvlTree tree = new AvlTree();

        // These variables will be modified by building and testing the AVL tree.
        int treeHeight = 0;
        long numberOfNodes = 0;

        // Build and test the AVL tree for the specified number of iterations.
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

            // Insert each coordinate into the AVL tree.
            long iTime = System.currentTimeMillis();
            for (int i = 0; i < coordinates.length; ++i) {
                if (tree.insert(coordinates[i])) {
                    // Create a KdNode that contains a TreeSet.
                    if (tree.insertedNode != null) {
                        final KdNode kdTreeNode = new KdNode(coordinates[i]);
                        tree.insertedNode.kdTreeNode = kdTreeNode;
                        kdTreeNode.avlTreeNode = tree.insertedNode;
                    }
                    // Verify the AVL tree after each insertion if requested.
                    if (verify) {
                        tree.verifyAvlTree();
                    }
                } else {
                    throw new RuntimeException("\n\nfailed to insert tuple " + i + "\n");
                }
            }
            iTime = System.currentTimeMillis() - iTime;
            insertTime[k] += (double) iTime / Constants.MILLISECONDS_TO_SECONDS;

            // Verify correct order of each node in the k-d tree and count the nodes.
            long vTime = System.currentTimeMillis();
            numberOfNodes = tree.verifyAvlTree();
            vTime = System.currentTimeMillis() - vTime;
            verifyTime[k] += (double) vTime / Constants.MILLISECONDS_TO_SECONDS;

            if (numberOfNodes + extraPoints != coordinates.length) {
                throw new RuntimeException("\n\nnumber of coordinates = " + coordinates.length + "  !=  " +
                                           "number of nodes + extra points = " + (numberOfNodes + extraPoints) + "\n");
            } else {
                treeHeight = tree.getHeight();
            }

            // Search for each coordinate in the AVL tree. 
            long sTime = System.currentTimeMillis();
            for (int i = 0; i < coordinates.length; ++i) {
                if (!tree.contains(coordinates[i])) {
                    throw new RuntimeException("\n\nfailed to find tuple " + i + "\n");
                }
            }
            sTime = System.currentTimeMillis() - sTime;
            searchTime[k] += (double) sTime / Constants.MILLISECONDS_TO_SECONDS;

            // Erase each coordinate from the AVL tree, and reverse
            // the order of the coordinates if reverse is true.
            // Don't reshuffle the coordinates so that the effect
            // of in-order and reverse-order erasure can be observed.
            long eTime = System.currentTimeMillis();
            if (reverse) {
                for (int i = coordinates.length - 1; i >= 0; --i) {
                    if (tree.erase(coordinates[i])) {
                        // Verify the AVL tree after each erasure if requested.
                        if (verify)
                        {
                            tree.verifyAvlTree();
                        }
                        // A search for a tuple after erasing it should fail. But note that erasure could
                        // remove a value from the TreeSet instead of removing the tuple from the AVL tree.
                        if (find && tree.contains(coordinates[i]) &&
                              tree.find(coordinates[i].getKey()).kdTreeNode.values.contains(coordinates[i].getValue())
                           )
                        {
                            throw new RuntimeException("\n\nfound tuple after erasing tuple " + i + "\n");
                        }
                        // A search for the prior tuple after erasing a tuple should succeed.
                        if (find && i > 0 && !tree.contains(coordinates[i-1]))
                        {
                            throw new RuntimeException("\n\nfailed to find prior tuple after erasing tuple " + i + "\n");
                        }
                    } else {
                        throw new RuntimeException("\n\nfailed to erase tuple " + i + "\n");
                    }
                }
            } else {
                for (int i = 0; i < coordinates.length; ++i) {
                    if (tree.erase(coordinates[i])) {
                        // Verify the AVL tree after each erasure if requested.
                        if (verify)
                        {
                            tree.verifyAvlTree();
                        }
                        // A search for a tuple after erasing it should fail. But note that erasure could
                        // remove a value from the TreeSet instead of removing the tuple from the AVL tree.
                        if (find && tree.contains(coordinates[i]) &&
                              tree.find(coordinates[i].getKey()).kdTreeNode.values.contains(coordinates[i].getValue())
                           )
                        {
                            throw new RuntimeException("\n\nfound tuple after erasing tuple " + i + "\n");
                        }
                        // A search for the next tuple after erasing a tuple should succeed.
                        if (find && i < coordinates.length-1 && !tree.contains(coordinates[i+1]))
                        {
                            throw new RuntimeException("\n\nfailed to find next tuple after erasing tuple " + i + "\n");
                        }
                    } else {
                        throw new RuntimeException("\n\nfailed to erase tuple " + i + " in dynamic tree\n");
                    }
                }
            }
            eTime = System.currentTimeMillis() - eTime;
            eraseTime[k] += (double) eTime / Constants.MILLISECONDS_TO_SECONDS;

            if ( !tree.isEmpty() ) {
                throw new RuntimeException("\n\ntree is not empty\n");
            }

            System.out.println("finished iteration " + (k + 1));
        }

        // Report AVLtree statistics.
        System.out.println("\nnumber of nodes = " + numberOfNodes + "  k-d tree height = " +
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

        // Report the rotation counters.
        System.out.println("\ninsert\tLL = " + (tree.lli/iterations) + "\tLR = " + (tree.lri/iterations) +
                           "\tRL = " + (tree.rli/iterations) + "\tRR = " + (tree.rri/iterations));
        System.out.println("delete\tLL = " + (tree.lle/iterations) + "\tLR = " + (tree.lre/iterations) +
                           "\tRL = " + (tree.rle/iterations) + "\tRR = " + (tree.rre/iterations));
    
        System.out.println();
    }
} // class TestKdTreeDynamic
