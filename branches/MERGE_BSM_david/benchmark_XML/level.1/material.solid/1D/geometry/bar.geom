*version
1.0
*title
1D bar/rod
*dimensions
6   # number of nodes
1   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   5   2
2   # number of node sets
# [ID] [nnd]
1   1
2   1
0   # number of side sets
*nodesets
*set
1   # number of nodes
1
*set
1   # number of nodes
6
# end node sets
*sidesets
*elements
*set
5   # number of elements
2   # number of element nodes
1  1  2  
2  2  3
3  3  4
4  4  5
5  5  6
# end elements
*nodes
6  # number of nodes
1  # number of spatial dimensions
 1  0.0
 2  1.0
 3  2.0
 4  3.0
 5  4.0
 6  5.0
 