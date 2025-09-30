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
import java.util.LinkedList;

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The NearestNeighborList class implements a fixed length list of both containing both a KdNode and euclidean distance
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
public class NearestNeighborList {
    protected long[] query; // point for the which the nearest neighbors will be found
    private int reqDepth; // requested number of nearest neighbors and therefore size of the above arrays
    protected KdNode[] nodes;  // set of nodes that are the nearest neighbors
    protected BigInteger[] dists; // set of distances from the query point to the above nodes
    protected int curDepth; // number of nearest nodes/distances on the list
    protected BigInteger curMaxDist; // distance to the last (and farthest) KdNode on the list
    protected boolean enable[];  // enable or disable each of the k dimensions

    /**
     * <p>
     * The {@code NearestNeighborList} constructor for NearestNeighborList class
     * </p>
     *
     * @param query - array containing the point to which the nearest neighbors will be found
     * @param numNeighbors - number of nearest neighbors to be found and hence size of the list
     */
    public NearestNeighborList( final long query[], final int numNeighbors) {
        this.nodes = new KdNode[numNeighbors];  // list of KdNodes
        this.dists = new BigInteger[numNeighbors];  // corresponding list of distances
        this.reqDepth = numNeighbors;
        this.curDepth = 0;  // the list is initially empty
        this.curMaxDist = BigInteger.valueOf(Long.MAX_VALUE).multiply(BigInteger.valueOf(Long.MAX_VALUE));
        this.query = new long[query.length];
        this.enable = new boolean[query.length];
       for (int i = 0; i < query.length; i++) {
            this.query[i] = query[i];
            this.enable[i] = true;
        }
    }

    public NearestNeighborList(final long query[], final int numNeighbors, final boolean enable[]) {
        this.nodes = new KdNode[numNeighbors];  // list of KdNodes
        this.dists = new BigInteger[numNeighbors];  // corresponding list of distances
        this.reqDepth = numNeighbors;
        this.curDepth = 0;  // the list is initially empty
        this.curMaxDist = BigInteger.valueOf(Long.MAX_VALUE).multiply(BigInteger.valueOf(Long.MAX_VALUE));
        this.query = new long[query.length];
        this.enable = new boolean[query.length];
        for (int i = 0; i < query.length; i++) {
            this.query[i] = query[i];
            this.enable[i] = enable[i];
        }
    }

    /**
     * <p>
     * The {@code add} Adds a newNode to this NearestNeighborList if its distance to the
     * query point is less than curMaxDistance
     * </p>
     *
     * @param newNode - KdNode to potentially be added to the list
     */
    protected void add(final KdNode newNode) {

        // Calculate the distance squared by subtracting the query coordinate
        // from the tuple coordinate and summing the squared differences.
        BigInteger dist = BigInteger.ZERO;
        for (int i = 0; i < newNode.tuple.length; i++) {
            if (enable[i]) {
                final long coor = newNode.tuple[i] - query[i];
                final BigInteger square =
                    BigInteger.valueOf(coor).multiply(BigInteger.valueOf(coor));
                dist = dist.add(square);
            }
        }
        // Check to see if this tuple should be added to the list.
        if ( dist.compareTo(curMaxDist) < 0 ) {
            // Need to add the tuple to the list, so do an insertion sort;
            // first do a binary search to find the index (pos) where the
            // new value is to be inserted.
            int pos = 0;
            if (curDepth > 0 && dist.compareTo(dists[0]) > 0) {
                int stride = (Integer.highestOneBit(curDepth-1) << 1);
                for (; stride > 0 && curDepth > 0; stride >>= 1) {
                    int newPos = pos + stride;
                    newPos = newPos > curDepth - 1 ? curDepth - 1 : newPos;
                    if (dist.compareTo(dists[newPos]) > 0) {
                        pos = newPos;
                    }
                }
                pos++;
            }
            // Next move the stuff in and after the insertion point down one place...
            if (curDepth < reqDepth) curDepth++;
            for (int i = curDepth-1;  i > pos; i--) {
                nodes[i] = nodes[i-1];
                dists[i] = dists[i-1];
            }
            // ...and put the new data in the open spot
            nodes[pos] = newNode;
            dists[pos] = dist;
            // if the list is full, replace the current max distance of the list to the last entry.
            if (curDepth == reqDepth) {
                curMaxDist = dists[reqDepth-1];
            }
        }
        return;
    }
} // class NearestNeighborList
