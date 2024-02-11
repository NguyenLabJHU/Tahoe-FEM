*version
1.0
*title
3 x 3 element square patch
*dimensions
16  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   9   4
4   # number of node sets
# [ID] [nnd]
1   4
2   4
3   1
4   1
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  3  4
*set
4   # number of nodes
13 14 15 16
*set
1   # number of nodes
1
*set
1   # number of nodes
13
# end node sets
*sidesets
*elements
*set
9  # number of elements
4   # number of element nodes
       1   1   2   6   5
       2   2   3   7   6
       3   3   4   8   7
       4   5   6  10   9
       5   6   7  11  10
       6   7   8  12  11
       7   9  10  14  13
       8  10  11  15  14 
       9  11  12  16  15
# end elements
*nodes
16  # number of nodes
2   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   1.5000000e+00   0.0000000e+00
    3   3.0000000e+00   0.0000000e+00
    4   4.5000000e+00   0.0000000e+00
    5   0.0000000e+00   1.5000000e+00
    6   1.5000000e+00   1.5000000e+00
    7   3.0000000e+00   1.5000000e+00
    8   4.5000000e+00   1.5000000e+00
    9   0.0000000e+00   3.0000000e+00
   10   1.5000000e+00   3.0000000e+00
   11   3.0000000e+00   3.0000000e+00
   12   4.5000000e+00   3.0000000e+00
   13   0.0000000e+00   4.5000000e+00
   14   1.5000000e+00   4.5000000e+00
   15   3.0000000e+00   4.5000000e+00
   16   4.5000000e+00   4.5000000e+00
