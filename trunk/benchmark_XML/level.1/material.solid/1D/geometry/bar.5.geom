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
6   # number of node sets
# [ID] [nnd]
1   1
2   1
3   1
4   1
5   1
6   1
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
*set
1   # number of nodes
5
*set
1   # number of nodes
6
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
 2  0.0
 3  1.0
 4  2.0
 5  3.0
 6  4.0
 1  5.0
 
