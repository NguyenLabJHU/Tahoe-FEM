*version
1.0
*title
2 linear quads with CSE between
*dimensions
10  # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   5
2   1   6
4   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   1 # LL
4   1 # UL
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
1   # number of nodes
1
*set
1   # number of nodes
7
# end node sets
*sidesets
*elements
*set
2   # number of elements
5   # number of element nodes
1  4  3  1  2  9
2  5  6  8  7  10
*set
1   # number of elements
6   # number of element nodes
1  3  4  6  5  9  10
# end elements
*nodes
10 # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  0.0  1.0
4  1.0  1.0
5  0.0  1.0
6  1.0  1.0
7  0.0  2.0
8  1.0  2.0
9  0.5  1.0
10 0.5  1.0