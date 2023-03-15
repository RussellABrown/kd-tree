The k-d tree-building algorithms are described in the following online articles.

http://www.jcgt.org/published/0004/01/03/

https://arxiv.org/abs/1410.5420

The Journal of Computer Graphics Techniques (JCGT) article contains a detailed description of an O(kn log n) k-d tree building algorithm and compares the performance of that algorithm to the performance of an O(n log n) algorithm.

In addition to the description of the O(kn log n) algorithm provided by the JCGT article, the arXiv article includes an appendix that describes improvements to the O(kn log n) and O(n log n) algorithms that were implemented following the presentation of those algorithms in the Journal of Computer Graphics Techniques article.

The kdTreeKnlogn.cpp and kdTreeNlogn.cpp source-code files build a balanced k-d tree.

The kdTreeMapKnlogn.cpp, kdTreeMapNlogn.cpp, kdTreeKmapKnlog.cpp, and kdTreeKmapNlogn.cpp source-code files build a balanced k-d tree-based key-to-value map. The latter two implementations of the key-to-value map execute twice as rapidly as the former two implementations. See the arXiv article for details.

All six source-code files include algorithms that search a k-d tree (1) for all points that lie inside a k-dimensional hyper-rectangular region; (2) for the m nearest neighbors to a query point sorted according to their distances to the query point via a priority queue; and (3) for the reverse nearest neighbors to each point in the k-d tree, where the reverse nearest neighbors to a given point are defined as the set of points to which that point is a nearest neighbor.

All six source-code files build a k-d tree and search a k-dimensional hyper-rectangular region using multiple threads.

The kdTreeKnlogn.cpp and kdTreeNlogn.cpp source-code files search for nearest neighbors and reverse nearest neighbors via a single thread.

The kdTreeMapKnlogn.cpp, kdTreeMapNlogn.cpp, kdTreeKmapKnlog.cpp, and kdTreeKmapNlogn.cpp source-code files search for the nearest neighbors to a single point via a single thread but search for the nearest neighbors and reverse nearest neighbors to all points in the k-d tree via multiple threads.

The kdTreeKnlogn.cpp and kdTreeNlogn.cpp source-code files store the nearest-neighbor and reverse-nearest-neighbor lists as maps wherein the key is a k-d node pointer and the value is a neighbor list.

The kdTreeMapKnlogn.cpp, kdTreeMapNlogn.cpp, kdTreeKmapKnlog.cpp, and kdTreeKmapNlogn.cpp source-code files store the nearest-neighbor and reverse-nearest-neighbor lists as vectors wherein each element is a neighbor list. The vectors permit faster access than maps permit.

The command-line options for controlling execution of the k-d tree-building programs are as follows.

-n The number of randomly generated points used to build the k-d tree

-m The number of nearest neighbors kept on the priority queue created when searching the k-d tree for nearest neighbors

-x The number of duplicate points added to test removal of duplicate points

-d The number of dimensions (k) of the k-d tree; not applicable to the kdTreeKmapKnlogn.cpp and kdTreeKmapNlogn.cpp programs because for these programs the number of dimensions is specified at compile time

-t The number of threads used to build and region-search the k-d tree

-s The search distance used for region search

-p The maximum number of nodes to print when reporting k-d tree search results

-b Enable nearest neighbors by exhaustive search for comparison to k-d tree search

-c Enable region search by exhaustive search for comparison to k-d tree search

-r Enable construction of nearest-neighbors and reverse-nearest-neighbors maps
