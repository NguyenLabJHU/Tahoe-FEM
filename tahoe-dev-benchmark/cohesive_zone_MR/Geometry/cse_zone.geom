*version
1.0
*title
2 linear quads without CSE between
*dimensions
6   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   4
6   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   1 # LL
4   1 # UL
5   3 # Left Side
6   3 # Right Side
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
5  6
*set
1   # number of nodes
1
*set
1   # number of nodes
5
*set
3   # number of nodes
1  3  5 
*set
3   # number of nodes
2  4  6
# end node sets
*sidesets
# end side sets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  4  3
2  3  4  6  5
# end elements
*nodes
6  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  0.0  1.0
4  1.0  1.0
5  0.0  2.0
6  1.0  2.0

