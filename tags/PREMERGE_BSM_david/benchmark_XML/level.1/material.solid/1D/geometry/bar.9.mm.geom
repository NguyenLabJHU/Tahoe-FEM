*version
1.0
*title
1D bar/rod
*dimensions
10  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   8   2
2   1   2
10  # number of node sets
# [ID] [nnd]
1   1
2   1
3   1
4   1
5   1
6   1
7   1
8   1
9   1
10  1
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
*set
1   # number of nodes
7
*set
1   # number of nodes
8
*set
1   # number of nodes
9
*set
1   # number of nodes
10
# end node sets
*sidesets
*elements
*set
8   # number of elements
2   # number of element nodes
1  2  3
2  3  4
3  4  5
4  5  6
5  7  8
6  8  9
7  9  10
8  10 1
*set
1   # number of elements
2   # number of element nodes
1  6  7
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
