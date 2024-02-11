*version
1.0
*title
1D bar/rod
*dimensions
6   # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   4   2
2   1   2
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
4   # number of elements
2   # number of element nodes
1  2  3  
2  3  4
3  5  6
4  6  1
*set
1   # number of elements
2   # number of element nodes
1  4  5
# end elements
*nodes
6  # number of nodes
1  # number of spatial dimensions
 2    0.00000000
 3   20.00000000
 4   40.00000000
 5   60.00000000
 6   80.00000000
 1  100.00000000 
