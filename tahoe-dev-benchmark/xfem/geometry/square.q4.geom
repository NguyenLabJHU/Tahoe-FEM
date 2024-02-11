*version
1.0
*title
1 element square with 4 nodes
*dimensions
4  # number of nodes
2  # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   1   4

4   # number of node sets
# [ID] [nnd]
1   2
2   2
3   2
4   4

0   # number of side sets
# end dimensions

*nodesets
*set
2   # number of nodes
1 2
*set
2   # number of nodes
1 4
*set
2   # number of nodes
3 4
*set
4   # number of nodes
1 2 3 4
# end node sets

*sidesets

*elements
*set
1 # number of elements
4 # number of element nodes
1 1 2 3 4
# end elements

*nodes
4  # number of nodes
2  # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00
    3   1.0000000e+00   1.0000000e+00
    4   0.0000000e+00   1.0000000e+00
    