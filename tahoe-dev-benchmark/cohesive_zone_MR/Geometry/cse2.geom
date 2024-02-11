*version
1.0
*title
2 linear quads with CSE between
*dimensions
8   # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   4
2   1   4
6   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   1 # LL
4   1 # UL
5   2 # followers
6   2 # leaders
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  3
*set
2   # number of nodes
6  8
*set
1   # number of nodes
1
*set
1   # number of nodes
6
*set
2   # number of nodes
2  4 
*set
2   # number of nodes
5  7
# end node sets
*sidesets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  4  3
2  5  6  8  7
*set
1   # number of elements
4   # number of element nodes
1  4  2  5  7
# end elements
*nodes
8  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  0.0  1.0
4  1.0  1.0
5  1.0  0.0
6  2.0  0.0
7  1.0  1.0
8  2.0  1.0

