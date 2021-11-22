This Github repository comprises the source files “kdTreeKnlogn.cpp”, “kdTreeNlogn.cpp”,
"kdTreeMapKnlogn.cpp", "kdTreeMapNlogn.cpp", "kdTreeKmapKnlog.cpp", and "kdTreeKmapNlogn.cpp"
that are C++ implementations of the k-d tree-building algorithms described in the online
articles

http://www.jcgt.org/published/0004/01/03/

https://arxiv.org/abs/1410.5420

The arXiv article includes an appendix that describes improvements to the k-d
tree-building algorithms relative to the description published in the Journal of
Computer Graphics Techniques article in 2015.

The source files “kdTreeKnlogn.cpp” and “kdTreeNlogn.cpp” are implementations that
build a k-d tree.

The source files "kdTreeMapKnlogn.cpp", "kdTreeMapNlogn.cpp", "kdTreeKmapKnlog.cpp",
and "kdTreeKmapNlogn.cpp" are implementations that build a k-d tree-based key-to-value
map. These latter two implementation execute twice as rapidly as the former two
implementations. See the arXiv article for further details.

All six source files include algorithms that (1) search a k-d tree for all points
that lie inside a k-dimensional hyper-rectangular region and (2) search a k-d tree
for the n nearest neighbors to a query point and sort those nearest neighbors
according to their distances to the query point via a priority queue.
