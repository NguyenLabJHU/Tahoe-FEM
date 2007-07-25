*version
1.0
*title
1D bar/rod
*dimensions
3   # number of nodes
1   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   2
3   # number of node sets
# [ID] [nnd]
1   1
2   1
3   1
0   # number of side sets
*nodesets
*set
1   # number of nodes
2
*set
1   # number of nodes
1
*set
1   # number of nodes
3
# end node sets
*sidesets
*elements
*set
2   # number of elements
2   # number of element nodes
1  2  3
2  3  1
# end elements
*nodes
3   # number of nodes
1   # number of spatial dimensions
 2  0.000000000
 3  2.500000000
 1  5.000000000
