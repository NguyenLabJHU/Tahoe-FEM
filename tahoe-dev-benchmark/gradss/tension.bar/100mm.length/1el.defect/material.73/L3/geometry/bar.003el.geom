*version
1.0
*title
1D bar/rod
*dimensions
7   # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   3
2   1   3
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
1 2 4 3
2 6 1 7
*set
1   # number of elements
3   # number of element nodes
1 4 6 5
# end elements
*nodes
7   # number of nodes
1   # number of spatial dimensions
2	0
3	16.66666667
4	33.33333333
5	50
6	66.66666667
7	83.33333333
1	100
