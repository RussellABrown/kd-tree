This Github repository comprises source code that implements two algorithms
to build a balanced k-d tree and a k-d tree-based
key-to-value map.

Specifically, the source-code files  “kdTreeKnlogn.cpp”, “kdTreeNlogn.cpp”,
"kdTreeMapKnlogn.cpp", "kdTreeMapNlogn.cpp", "kdTreeKmapKnlog.cpp", and
"kdTreeKmapNlogn.cpp" are C++ implementations of the k-d tree-building
algorithms described in the following online articles.

http://www.jcgt.org/published/0004/01/03/

https://arxiv.org/abs/1410.5420

The Journal of Computer Graphics Techniques article contains a detailed
description of an O(kn logn) k-d tree building algorithm and compares the
performance of that algorithm to the performance of an O(n logn) algorithm.

The arXiv article includes an appendix that describes improvements to two
k-d tree-building algorithms that were implemented following the description
of those algorithms in the Journal of Computer Graphics Techniques article.

The source files “kdTreeKnlogn.cpp” and “kdTreeNlogn.cpp” are implementations
that build a balanced k-d tree.

The source files "kdTreeMapKnlogn.cpp", "kdTreeMapNlogn.cpp", "kdTreeKmapKnlog.cpp",
and "kdTreeKmapNlogn.cpp" are implementations that build a balanced k-d tree-based
key-to-value map. The latter two implementations of the key-to-value map execute
twice as rapidly as the former two implementations. See the arXiv article for details.

All six source files include algorithms that (1) search a k-d tree for all points
that lie inside a k-dimensional hyper-rectangular region and (2) search a k-d tree
for the n nearest neighbors to a query point and sort those n nearest neighbors
by increasing distance to the query point via a priority queue.

The command-line options for controlling execution of the k-d tree-building
programs are as follows.

-n The number of randomly generated points used to build the k-d tree

-m The number of nearest neighbors kept on the priority queue created
   when searching the k-d tree for nearest neighbors
   
-x The number of duplicate points added to test removal of duplicate points

-d The number of dimensions (k) of the k-d tree; not applicable to the
   kdTreeKmapKnlogn.cpp and kdTreeKmapNlogn.cpp programs because for these
   programs the number of dimensions is specified at compile time
   
-t The number of threads used to build and region-search the k-d tree

-s The search distance used for region search

-p The maximum number of nodes to print when reporting k-d tree search results

-b Enable nearest neighbors by exhaustive search for comparison to k-d tree search

-c Enable region search by exhaustive search for comparison to k-d tree search

