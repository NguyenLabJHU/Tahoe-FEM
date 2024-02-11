*version
1.0
*title
1D bar/rod
*dimensions
7   # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   4   2
2   2   2
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
3  6  7
4  7  1
*set
2   # number of elements
2   # number of element nodes
1  4  5
2  5  6
# end elements
*nodes
7   # number of nodes
1   # number of spatial dimensions
2 0
3 16.67
4 33.33
5 50
6 66.67
7 83.33
1 100
