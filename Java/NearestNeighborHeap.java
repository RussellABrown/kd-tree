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

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The NearestNeighborHeap class implements a fixed length list of both containing both a KdNode and euclidean distance
 * from the tuple in the node to a query point.  When a KdNode is added to the list it is unconditionally placed in
 * the list until the list is full.  After the list is full, a KdNode is added to the list only if the calculated
 * distance from the query point to the tuple is less than the farthest KdNode currently in the list; and in that
 * case, the current farthest KdNode and distance are removed from the list to make room for it.
 *
 * The list is maintained in two corresponding fixed length arrays, one for the KdNodes and one for the distance to
 * the query point.  These arrays are stored in order of increasing distance.  When a new node is added, regardless
 * of whether or not the list is full, an insertion sort is done to place the new KdNode in the proper order in the list.
 * A binary search is done to find the correct location, then the farther entries are moved down to make room for the new
 * one, which results in the old farthest KdNode being dropped.  In this way, the farthest KdNode is always at the end of
 * the list when the list is full.
 *
 * A separate variable, curMaxDist, holds the value of the distance to compar against to decide whether to accept
 * a new KdNode.  This variable is initialized to a maximum omt value in the constructor so that all KdNodes will
 * be added to the list until the list is full.  Once the list is full, curMaxDist it set to the distance of the
 * last KdNode on the list, which is the farthest KdNode on the list, to enable distance comparison to that KdNode.
 * <p>
 */
public class NearestNeighborHeap {
    protected long[] query; // query point for which nearest neighbors will be found
    protected boolean[] enable;
    private int reqDepth; // requested number of nearest neighbors
    protected KdNode[] nodes; // array of pointers to KdNodes that are the nearest neighbors
    protected BigInteger[] dists; // vector of squared distances
    protected int curDepth; // number of nearest nodes/distances on the heap

    /*
    * Constructor that enables distance test on all dimensions
    *
    * Calling parameters:
    *
    * query - a vector that defines the query point
    * numNeighbors - the number of nearest neighbors desired
    */

    public NearestNeighborHeap(final long[] query, final int numNeighbors) {
        this.nodes = new KdNode[numNeighbors + 1]; // heap of KdNode* (address 0 is unused)
        this.dists = new BigInteger[numNeighbors + 1]; // corresponding heap of distances
        for (int i = 0; i < dists.length; ++i) {
            dists[i] = BigInteger.ZERO; // initialize each distance to 0
        }
        this.reqDepth = numNeighbors;
        this.curDepth = 0;
        this.query = new long[query.length];
        this.enable = new boolean[query.length];
        for (int i = 0; i < query.length; i++) {
            this.query[i] = query[i];
            this.enable[i] = true;
        }
    }

    /*
    * Constructor that enables distance test for only specified dimensions
    *
    * Calling parameters:
    *
    * query - a vector that defines the query point
    * numNeighbors - the number of nearest neighbors desired
    * enable - a vector that specifies the dimensions for which to test distance
    */

    public NearestNeighborHeap(final long[] query, final int numNeighbors, final boolean[] enable) {
        this.nodes = new KdNode[numNeighbors + 1]; // heap of KdNode* (address 0 is unused)
        this.dists = new BigInteger[numNeighbors + 1]; // corresponding heap of distances
        for (int i = 0; i < dists.length; ++i) {
            dists[i] = BigInteger.ZERO; // initialize each distance to 0
        }
        this.reqDepth = numNeighbors;
        this.curDepth = 0;
        this.query = new long[query.length];
        this.enable = new boolean[query.length];
        for (int i = 0; i < query.length; i++) {
            this.query[i] = query[i];
            this.enable[i] = enable[i];
        }
    }

    /*
    * Swap two elements in the heap.
    *
    * Calling parameters:
    *
    * i - the index of the first element
    * j - the index of the second element
    */
    private void swap(final int i, final int j) {
        
        final BigInteger tempDist = dists[i];
        final KdNode tempNode = nodes[i];
        dists[i] = dists[j];
        nodes[i] = nodes[j];
        dists[j] = tempDist;
        nodes[j] = tempNode;
    }

    /*
    * Allow an element to rise upward through the heap.
    *
    * Calling parameter:
    *
    * kk - the index of the element
    */
    private void rise(final int kk) {

        int k = kk;
        while (k > 1 && dists[k/2].compareTo(dists[k]) < 0) {
            swap(k/2, k);
            k = k/2;
        }
    }

    /*
    * Allow an element to fall downward through the heap.
    *
    * Calling parameter:
    *
    * kk - the index of the element
    */
    private void fall(final int kk) {

        int k = kk;
        while (2*k <= curDepth) {
            int j = 2*k;
            if (j < curDepth && dists[j].compareTo(dists[j+1]) < 0) {
                ++j;
            }
            if (dists[k].compareTo(dists[j]) >= 0) {
                break;
            }
            swap(k, j);
            k = j;
        }
    }



    /*
    * Remove the top element of the heap and re-order the remaining elements.
    *
    * return a pair that contains a pointer to the top KdNode and the distance to that KdNode
    */
    protected Paire removeTop() {
        Paire returnPair = new Paire(dists[1], nodes[1]);
        swap(1, curDepth--);
        nodes[curDepth+1] = null;
        fall(1);
        return returnPair;
    }

    /*
    * Add a new KdNode to the NearestNeighborHeap if its distance to the
    * query point is less than curMaxDistance.
    *
    * Calling parameter:
    *
    * node - KdNode to potentially add to the heap
    */
    protected void add(final KdNode node) {
        // Find the distance by subtracting the query from the tuple and
        // calculating the sum of the squared distances. Note that conversion
        // from type K to cpp_int may result in loss of precision but avoids
        // the possibility of integer overflow.
        BigInteger dist = BigInteger.ZERO;
        for (int i = 0; i < node.tuple.length; ++i) {
            if (enable[i]) {
                final BigInteger tup = BigInteger.valueOf(node.tuple[i]);
                final BigInteger que = BigInteger.valueOf(query[i]);
                final BigInteger coor = tup.subtract(que);
                final BigInteger square = coor.multiply(coor);
                dist = dist.add(square);
            }
        }
        // If the queue is not full, add the point to the bottom of the heap unconditionally and let it rise;
        if (!heapFull()) {
            dists[++curDepth] = dist;
            nodes[curDepth] = node;
            rise(curDepth);
        }
        // otherwise, if the point is closer than the top of the heap, overwrite the top and let it fall.
        else if (dist.compareTo(curMaxDist()) < 0) {
            dists[1] = dist;
            nodes[1] = node;
            fall(1);
        }
        return;
    }

    /* Return the current maximum distance, i.e., dists[1] */
    protected BigInteger curMaxDist() {
        return dists[1];
    }

    /* Return true if the heap is full. */
    protected boolean heapFull() {
        return curDepth >= reqDepth;
    }

    /* Return the current depth of the heap, i.e., the number of nearest nodes/distances elements on the heap. */
    protected int heapDepth() {
        return curDepth;
    }

} // class NearestNeighborHeap
