*version
1.0
*title
1 x 1 x 1 element patch
*dimensions
8   # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   8
4   # number of node sets
# [ID] [nnd]
1   8
2   4
3   4
4   1
0   # number of side sets
# end dimensions
*nodesets
*set
8   # number of nodes
2  3  5  4  
6  7  1  8
*set
4   # number of nodes
4  5  8  1
*set
4   # number of nodes
2  3  6  7
*set
1   # number of nodes
2
# end node sets
*sidesets
*elements
*set
1  # number of elements
8   # number of element nodes
1  2  3  5  4  6  7  1  8
# end elements
*nodes
8   # number of nodes
3   # number of spatial dimensions
    2   0.0000000e+00   0.0000000e+00   0.0000000e+00
    3   5.0000000e+01   0.0000000e+00   0.0000000e+00
    4   0.0000000e+00   1.0000000e+02   0.0000000e+00
    5   5.0000000e+01   1.0000000e+02   0.0000000e+00

    6   0.0000000e+00   0.0000000e+00   1.0000000e+01
    7   5.0000000e+01   0.0000000e+00   1.0000000e+01
    8   0.0000000e+00   1.0000000e+02   1.0000000e+01
    1   5.0000000e+01   1.0000000e+02   1.0000000e+01
