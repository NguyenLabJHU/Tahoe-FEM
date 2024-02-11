*version
1.0
*title
1x2 elements
*dimensions
6  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   4
8   # number of node sets
# [ID] [nnd]
1   2
2   3
3   2
4   3
5   1
6   1
7   1
8   1
0   # number of side sets
# end dimensions

*nodesets
*set
2   # number of nodes
1  4  
*set
3   # number of nodes
4  5  6  
*set
2   # number of nodes
3  6  
*set
3   # number of nodes
1  2  3  
*set
1   # number of nodes
1
*set
1   # number of nodes
4
*set
1   # number of nodes
6
*set
1   # number of nodes
3
# end node sets
*sidesets
*elements
*set
2  # number of elements
4   # number of element nodes
       1       1       4       5       2
       2       2       5       6       3
# end elements
*nodes
6 # number of nodes
2  # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   0.0000000e+00   1.0000000e+00
    3   0.0000000e+00   2.0000000e+00
    4   1.0000000e+00   0.0000000e+00
    5   1.0000000e+00   1.0000000e+00
    6   1.0000000e+00   2.0000000e+00
