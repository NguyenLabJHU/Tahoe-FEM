*version
1.0
*title
2 cubes stacked along the z-axes
*dimensions
12  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   8
2   # number of node sets
# [ID] [nnd]
1   4
2   4
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  3  4
*set
4   # number of nodes
9 10 11 12
# end node sets
*sidesets
*elements
*set
2  # number of elements
8  # number of element nodes
1  1  2  3  4  5  6  7  8
2  5  6  7  8  9 10 11 12
# end elements
*nodes
12  # number of nodes
3   # number of spatial dimensions
1  0.0  0.0  0.0
2  0.0  0.0  1.0
3  1.0  0.0  1.0
4  1.0  0.0  0.0

5  0.0  1.0  0.0
6  0.0  1.0  1.0
7  1.0  1.0  1.0
8  1.0  1.0  0.0

 9  0.0  2.0  0.0
10  0.0  2.0  1.0
11  1.0  2.0  1.0
12  1.0  2.0  0.0
