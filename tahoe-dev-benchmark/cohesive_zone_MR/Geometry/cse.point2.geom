*version
1.0
*title
1 CSE only
*dimensions
4   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   4
6   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   1 # LL
4   1 # UL
5   1 # LR
6   1 # UR
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
3  4
*set
1   # number of nodes
1
*set
1   # number of nodes
3
*set
1   # number of nodes
2
*set
1   # number of nodes
4
# end node sets
*sidesets
*elements
*set
1   # number of elements
4   # number of element nodes
1  1  2  4  3
# end elements
*nodes
4  # number of nodes
2  # number of spatial dimensions
1  0.0  1.0
2  0.4  0.0
3  0.0  1.0
4  0.4  0.0
