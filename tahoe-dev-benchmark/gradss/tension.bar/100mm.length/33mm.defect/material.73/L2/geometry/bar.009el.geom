*version
1.0
*title
1D bar/rod
*dimensions
10  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   6   2
2   3   2
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
6   # number of elements
2   # number of element nodes
1  2  3
2  3  4
3  4  5
4  8  9
5  9  10
6  10 1
*set
3   # number of elements
2   # number of element nodes
1  5  6
2  6  7
3  7  8
# end elements
*nodes
10  # number of nodes
1   # number of spatial dimensions
2    0.00000000
3   11.11111111
4   22.22222222
5   33.33333333
6   44.44444444
7   55.55555556
8   66.66666667
9   77.77777778
10  88.88888889
1  100.00000000
