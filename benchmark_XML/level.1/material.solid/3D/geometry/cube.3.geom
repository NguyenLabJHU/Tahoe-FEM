*version
1.0
*title
1 element cube
*dimensions
27  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   8
6   # number of node sets
# [ID] [nnd]
1   4
2   4
3   4
4   4
5   4
6   4
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  3  4
*set
4   # number of nodes
5  6  7  8
*set
4   # number of nodes
 1  2  5  6
*set
4   # number of nodes
 1  3  5  7
*set
4   # number of nodes
 2  4  6  8
*set
4   # number of nodes
 3  4  7  8
# end node sets
*sidesets
*elements
*set
cube.3.elem
# end elements
*nodes
cube.3.node
