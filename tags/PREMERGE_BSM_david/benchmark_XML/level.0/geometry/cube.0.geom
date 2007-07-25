*version
1.0
*title
1 cubic element
*dimensions
8   # number of nodes
3   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   1   8
2  16   2
3   # number of node sets
# [ID] [nnd]
1   8
2   4
3   1
0   # number of side sets
# end dimensions
*nodesets
*set
8   # number of nodes
1 2 3 4 5 6 7 8
*set
4
1 2 3 4
*set
1
7
# end node sets
*sidesets
*elements
*set
1   # number of elements
8   # number of element nodes
1  1  2  3  4  5  6  7  8
*set
16  # number of elements
2   # number of element nodes
1  1  2
2  2  3
3  3  4
4  4  1
5  5  6
6  6  7
7  7  8
8  8  5
9  1  5
10 2  6
11 3  7
12 4  8
13 1  7
14 2  8
15 3  5
16 4  6
# end elements
*nodes
8  # number of nodes
3  # number of spatial dimensions
 1  0.0  0.0  0.0
 2  1.0  0.0  0.0
 3  1.0  1.0  0.0
 4  0.0  1.0  0.0
 5  0.0  0.0  1.0
 6  1.0  0.0  1.0
 7  1.0  1.0  1.0
 8  0.0  1.0  1.0
