*version
1.0
*title
2 tetrahedron elements with a CSE between
*dimensions
8   # number of nodes
3   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   4
2   1   6
4   # number of node sets
# [ID] [nnd]
1   1 # bottom
2   1 # top
3   1 # fixed line_X
4   2 # fixed line_Y
0   # number of side sets
# end dimensions
*nodesets
*set
1   # number of nodes
8
*set
1   # number of nodes
4
*set
1   # number of nodes
5
*set
2   # number of nodes
5 6
# end node sets
*sidesets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  3  4
2  5  6  7  8
*set
1   # number of elements
6   # number of element nodes
1 5 6 7 1 2 3
# end elements
*nodes
8  # number of nodes
3  # number of spatial dimensions
1  0.00000  0.00000 0.00000
2  1.00000  0.00000 0.00000
3  0.50000  0.86603 0.00000
4  0.50000  0.57735 0.81650
5  0.00000  0.00000 0.00000
6  1.00000  0.00000 0.00000
7  0.50000  0.86603 0.00000
8  0.50000  0.57735 -0.81650

