The k-d tree-building algorithms are described in the following online articles.

http://www.jcgt.org/published/0004/01/03/

https://arxiv.org/abs/1410.5420

The Journal of Computer Graphics Techniques (JCGT) article contains a detailed description of an O(kn log n) k-d tree building algorithm and compares the performance of that algorithm to the performance of an O(n log n) algorithm.

In addition to the description of the O(kn log n) algorithm provided by the JCGT article, the arXiv article includes an appendix that describes improvements to the O(kn log n) and O(n log n) algorithms that were implemented following the publication of those algorithms in the JCGT article.

The test_kdtree.cpp file and associated .h files build and test a k-d tree, which implements a set, via either the O(kn log n) or the O(n log n) algorithm.

The test_kdmap.cpp file and associated .h files build and test a k-d tree-based key-to-multiple-value map via either the O(kn log n) or the O(n log n) algorithm.

The '-D NLOGN' compilation define specifies the O(n log n) algorithm; otherwise, the O(kn log n) algorithm is used.

The '-D PREALLOCATE' compilation define causes temporary data structures that are required to build the k-d tree to be allocated and deallocated en masse, which improves performance.

For test_kdmap.cpp, the '-D DIMENSIONS=k' compilation define specifies the number of dimensions k. This define is useful only if the '-D PREALLOCATE' define fails to compile correctly and is ignored unless '-D PREALLOCATE' is specified as well. It improves performance similarly to '-D PREALLOCATE' but requires that the number of dimensions be specified at compile time instead of run time, so it results in less flexibility than '-D PREALLOCATE'.

Also for for test_kdmap.cpp,the '-D TREE' compilation define causes a k-d tree to be built instead of a k-d tree-based key-to-multiple-value map.

See the test_kdtree.cpp and test_kdmap.cpp files for discussion of other compilation defines.

The command-line options that control execution of the main() function to build and search the k-d tree are as follows:

-i The number of iterations of k-d tree creation; more iterations enable more reliable measurement of execution times                                                    
-n The number of randomly generated points used to build the k-d tree

-m The number of nearest neighbors kept on the priority queue created when searching the k-d tree for nearest neighbors

-x The number of duplicate points added to test removal of duplicate points

-d The number of dimensions (k) of the k-d tree; not applicable to the test_kdset.cpp program if '-D DIMENSIONS=k' is specified to the compiler

-t The number of threads used to build and region-search the k-d tree

-s The search distance used for region search

-p The maximum number of nodes to print when reporting k-d tree search results

-b Enable nearest neighbors by exhaustive search for comparison to k-d tree search

-c Enable region search by exhaustive search for comparison to k-d tree search

-r Enable construction of nearest-neighbors and reverse-nearest-neighbors maps or vectors
