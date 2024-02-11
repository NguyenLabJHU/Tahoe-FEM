*version
1.0
*title
1D bar/rod
*dimensions
5   # number of nodes
1   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   3
2   # number of node sets
# [ID] [nnd]
1   1
2   1
0   # number of side sets
*nodesets
*set
1   # number of nodes
2
*set
1   # number of nodes
1
# end node sets
*sidesets
*elements
*set
2   # number of elements
3   # number of element nodes
1  2  4  3
2  4  1  5
# end elements
*nodes
5   # number of nodes
1   # number of spatial dimensions
 2  0.000000000
 3  25.00000000
 4  50.00000000
 5  75.00000000
 1  100.0000000
