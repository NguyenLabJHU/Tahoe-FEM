*version
1.0
*title
4 triangle elements with 4 CSE between
*dimensions
12  # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   4   3
2   4   4
6   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   2 # LL
4   2 # RL
5   1 # LBN
6   1 # RBN
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
7  8
*set
2   # number of nodes
10  11
*set
2   # number of nodes
4  5
*set
1   # number of nodes
1
*set
1   # number of nodes
2
# end node sets
*sidesets
*elements
*set
4   # number of elements
3   # number of element nodes
1  1  2  3
2  4  5  6
3  7  8  9
4  10  11  12
*set
4   # number of elements
4   # number of element nodes
1  3  2  4  6
2  6  5  7  9
3  9  8  10  12
4  12  11  1  3
# end elements
*nodes
12  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  0.4  0.7
4  1.0  0.0
5  1.0  1.0
6  0.4  0.7
7  1.0  1.0
8  0.0  1.0
9  0.4  0.7
10  0.0  1.0
11  0.0  0.0
12  0.4  0.7
