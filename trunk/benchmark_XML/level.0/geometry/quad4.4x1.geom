*version
1.0
*title
4 x 1 patch of quad 4 elements
*dimensions
10  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   4
4   # number of node sets
# [ID] [nnd]
1   5
2   5
3   2
4   2
1   # number of side sets
# [ID] [element set ID] [ns]
1   1   1
# end dimensions
*nodesets
*set
5   # number of nodes
1  2  3  4  5
*set
5   # number of nodes
6  7  8  9  10
*set
2   # number of nodes
1  6
*set
2   # number of nodes
5  10
# end node sets
*sidesets
*set
1   # number of sides
# [element] [face]
4  2 
# end side sets
*elements
*set
4  # number of elements
4  # number of element nodes
1  1  2  7  6
2  2  3  8  7
3  3  4  9  8
4  4  5 10  9
# end elements
*nodes
10  # number of nodes
2   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00
    3   2.0000000e+00   0.0000000e+00
    4   3.0000000e+00   0.0000000e+00
    5   4.0000000e+00   0.0000000e+00

    6   0.0000000e+00   1.0000000e+00
    7   1.0000000e+00   1.0000000e+00
    8   2.0000000e+00   1.0000000e+00
    9   3.0000000e+00   1.0000000e+00
   10   4.0000000e+00   1.0000000e+00
