*version
1.0
*title
1 x 1 x 2 element beam
*dimensions
12  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   8
4   # number of node sets
# [ID] [nnd]
1   4
2   1
3   1 # set with all model nodes 
4   8 # all nodes except those other node sets
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  3  4
*set
1   # number of nodes
11
*set
1
-1
*set
8
5 6 7 8 9 10 11 12
# end node sets
*sidesets
*elements
*set
2   # number of elements
8   # number of element nodes
1  1  2  3  4  5  6  7  8
2  5  6  7  8  9 10 11 12
# end elements
*nodes
12  # number of nodes
 3  # number of spatial dimensions
 1  0.0  0.0  0.0
 2  1.0  0.0  0.0
 3  1.0  1.0  0.0
 4  0.0  1.0  0.0
 5  0.0  0.0  1.0
 6  1.0  0.0  1.0
 7  1.0  1.0  1.0
 8  0.0  1.0  1.0
 9  0.0  0.0  2.0
10  1.0  0.0  2.0
11  1.0  1.0  2.0
12  0.0  1.0  2.0
