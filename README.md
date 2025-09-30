DYNAMIC K-D TREE:

The dynamic k-d tree-building algorithms are described in the following online article.

https://arxiv.org/abs/2509.08148

The dynamic k-d tree-building algorithms maintain a balanced k-d tree that resembles an AVL or red-black tree, and that self-balances by rebuilding a subtree that becomes unbalanced as a result of insertion or deletion of a k-dimensional tuple. The multi-threaded, static k-d tree-building algorithms described hereafter rebuild that subtree.

The dynamic k-d tree implements a set. Also included is the implementation of a dynamic k-d tree-based key-to-multiple-value map.

The dynamic k-d tree is implemented as the KdTreeDynamic class that is derived from the KdTree class that implements the static k-d tree. Hence, search algorithms implemented by the KdTree class, such as region search and nearest-neighbor search, may be performed for either a static k-d tree or a dynamic k-d tree. Also, it is possible to build a static k-d tree initially via one of the two the static algorithms, and thereafter modify that k-d tree via the dynamic algorithms.

The test_kdtreedynamic.cpp file and associated .h files (kdTreeDynamic.h, kdTreeNode.h, KdTreeNlogn.h, kdTreeKnlogn.h, kdTreeHeapSort.h and KdTreeMergeSort.h) build and test a dynamic k-d tree.

The test_kdtreedynamic.stats.cpp file and associated .h files kdTreeDynamic.stats.h, kdTreeNode.h, KdTreeNlogn.h, kdTreeKnlogn.h, kdTreeHeapSort.h and KdTreeMergeSort.h) build and test a dynamic k-d tree and optionally collect various performance statistics, including histograms.

The test_kdmapdynamic.cpp file and associated .h files (kdMapDynamic.h, kdMapNode.h, kdMapNlogn.h, kdMapKnlogn.h, kdMapHeapSort.h and kdMapMergeSort.h) build and test a dynamic k-d tree-based key-to-multiple-value map.

See the test_kdtreedynamic.cpp, test_kdtreedynamic.stats.cpp, and test_kdmapdynamic.cpp files for compilation instructions, explanation of -D compilation defines, and explanation of command-line options.

The TestKdTreeDynamic.java file and associated .java files (KdTreeDynamic.java, KdTree.java, KdTreeNlogn.java, KdTreeKnlogn.java, KdNode.java, MergeSort.java, NearestNeighborList.java, Pair.java and Constants.java) build and test a dynamic k-d tree-based key-to-multiple-value map. See the comments in TestKdTreeDynamic.java for an explanation of the configuration constants in Constants.java.

The KdTreeDynamic class extends the KdTree class.


STATIC K-D TREE:

The multi-threaded, static k-d tree-building algorithms are described in the following online articles.

http://www.jcgt.org/published/0004/01/03/

https://arxiv.org/abs/1410.5420

https://arxiv.org/abs/2506.20687

The Journal of Computer Graphics Techniques (JCGT) article contains a detailed description of an O[knlog(n)] k-d tree building algorithm and compares the performance of that algorithm to the performance of an O[nlog(n)] algorithm.

In addition to the description of the O[knlog(n)] algorithm provided by the JCGT article, the first of the two arXiv articles includes an appendix that describes improvements to the O[knlog(n)] and O[nlog(n)] algorithms that were implemented following publication of those algorithms in the JCGT article.

The second of the two arXiv articles is a review that contrasts the O[knlog(n)] and O[nlog(n)] algorithms to an algorithm proposed by Yu Cao et al. in the journal article, "A New Method to Construct the KD Tree Based on Presorted Results", Complexity Volume 2020, Article ID 8883945, https://doi.org/10.1155/2020/8883945

All three of the above-listed arXiv articles include C++ implementations of the improvements that followed the publication of those algorithms in the JCGT article and are described as follows. Those same improvements are available at this Github web page.

The test_kdtree.cpp file and associated .h files (kdTreeNode.h, kdTreeNlogn.h, kdTreeKnlogn.h, kdTreeYuCao.h, kdTreeHeapSort.h and kdTreeMergeSort.h) build and test a static k-d tree, which implements a set, via the O[knlog(n)], O[nlog(n)], or O[knlog(n)] + O[nlog(n)] algorithm.

The test_kdmap.cpp file and associated .h files (kdMapNode.h, kdMapNlogn.h, kdMapKnlogn.h, kdMapHeapSort.h and kdMapMergeSort.h) build and test a static k-d tree-based key-to-multiple-value map via either the O[knlog(n)] or the O[nlog(n)] algorithm.

See the test_kdtree.cpp and test_kdmap.cpp files for compilation instructions, explanation of -D compilation defines, and explanation of command-line options.

The TestKdTree.java file and associated .java files (KdTree.java, KdTreeNlogn.java, KdTreeKnlogn.java, KdNode.java, MergeSort.java, NearestNeighborList.java, Pair.java and Constants.java) build and test a static k-d tree-based key-to-multiple-value map. See the comments in TestKdTree.java for an explanation of the configuration constants in Constants.java.

