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
1
*set
1   # number of nodes
2
*set
1   # number of nodes
3
# end node sets
*sidesets
*elements
*set
2   # number of elements
2   # number of element nodes
1  1  2
2  2  3
# end elements
*nodes
3   # number of nodes
1   # number of spatial dimensions
 1    0.000000000
 2   50.000000000
 3  100.000000000
