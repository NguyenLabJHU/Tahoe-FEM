*version
1.0
*title
simple beam
*dimensions
12  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   5   4
2   # number of node sets
# [ID] [nnd]
1   2
2   1
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  7
*set
1   # number of nodes
12
# end node sets
*sidesets
*elements
*set
5  # number of elements
4  # number of element nodes
1  1  2  8  7
2  2  3  9  8
3  3  4  10  9
4  4  5  11  10
5  5  6  12  11
# end elements
*nodes
12  # number of nodes
2   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00
    3   2.0000000e+00   0.0000000e+00
    4   3.0000000e+00   0.0000000e+00
    5   4.0000000e+00   0.0000000e+00
    6   5.0000000e+00   0.0000000e+00

    7   0.0000000e+00   1.0000000e+00
    8   1.0000000e+00   1.0000000e+00
    9   2.0000000e+00   1.0000000e+00
   10   3.0000000e+00   1.0000000e+00
   11   4.0000000e+00   1.0000000e+00
   12   5.0000000e+00   1.0000000e+00
