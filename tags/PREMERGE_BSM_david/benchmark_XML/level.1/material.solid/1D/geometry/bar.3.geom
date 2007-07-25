*version
1.0
*title
1D bar/rod
*dimensions
4   # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   2
2   1   2
4   # number of node sets
# [ID] [nnd]
1   1
2   1
3   1
4   1
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
*set
1   # number of nodes
4
# end node sets
*sidesets
*elements
*set
2   # number of elements
2   # number of element nodes
1  2  3
2  4  1
*set
1   # number of elements
2   # number of element nodes
1  3  4
# end elements
*nodes
4   # number of nodes
1   # number of spatial dimensions
 2  0.000000000
 3  1.666666667
 4  3.333333333
 1  5.000000000
