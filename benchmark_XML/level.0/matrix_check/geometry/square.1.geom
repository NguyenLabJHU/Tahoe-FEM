*version
1.0
*title
2 x 2 element square patch
*dimensions
9   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   4
4   # number of node sets
# [ID] [nnd]
1   3
2   3
3   1
4   1
0   # number of side sets
# end dimensions
*nodesets
*set
3   # number of nodes
1  2  3
*set
3   # number of nodes
7  8  9
*set
1   # number of nodes
1
*set
1   # number of nodes
7
# end node sets
*sidesets
*elements
*set
4  # number of elements
4  # number of element nodes
1  1  2  5  4
2  2  3  6  5
3  4  5  8  7
4  5  6  9  8

# end elements
*nodes
9   # number of nodes
2   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00
    3   2.0000000e+00   0.0000000e+00

    4   0.0000000e+00   1.0000000e+00
    5   1.0000000e+00   1.0000000e+00
    6   2.0000000e+00   1.0000000e+00

    7   0.0000000e+00   2.0000000e+00
    8   1.0000000e+00   2.0000000e+00
    9   2.0000000e+00   2.0000000e+00
