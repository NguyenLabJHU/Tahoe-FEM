*version
1.0
*title
1 element cube 
*dimensions
8  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   8
7   # number of node sets
# [ID] [nnd]
1   4
2   4
3   4
4   4
5   4
6   4
7   2
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes (fixed in x direction)
1  3  5  7
*set
4   # number of nodes (fixed in y direction)
1  2  3  4
*set
4   # number of nodes (fixed in z direction)
3  4  7  8  
*set
4   # number of nodes (confining pressure in -x direction)
2  4  6  8
*set
4   # number of nodes (confining pressure in -y direction)
5  6  7  8
*set
4   # number of nodes (top: displacement /confining pressure in z direction)
1  2  5  6
*set
2   # number of nodes (top surface nodes free in -x direction)
2  6
# end node sets
*sidesets
*elements
*set
1  # number of elements
8   # number of element nodes
1  3  1  2  4  7  5  6  8
# end elements
*nodes
8  # number of nodes
3   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00    0.0000000e+00
    2   1.0000000e+00   0.0000000e+00    0.0000000e+00
    3   0.0000000e+00   0.0000000e+00   -1.0000000e+00
    4   1.0000000e+00   0.0000000e+00   -1.0000000e+00
    5   0.0000000e+00   1.0000000e+00    0.0000000e+00
    6   1.0000000e+00   1.0000000e+00    0.0000000e+00
    7   0.0000000e+00   1.0000000e+00   -1.0000000e+00
    8   1.0000000e+00   1.0000000e+00   -1.0000000e+00
   

