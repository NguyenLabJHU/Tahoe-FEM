*version
1.0
*title
2 x 2 patch of quad 8 elements
*dimensions
21  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   8
6   # number of node sets
# [ID] [nnd]
1   5
2   5
3   5
4   5
5   1
6   1
1   # number of side sets
# [ID] [element set ID] [ns]
1   1   2
# end dimensions
*nodesets
*set
5   # number of nodes
1  2  3  4  5
*set
5   # number of nodes
5  8 13 16 21
*set
5   # number of nodes
17 18 19 20 21
*set
5   # number of nodes
1  6  9 14 17
*set
1   # number of nodes
1
*set
1   # number of nodes
21
# end node sets
*sidesets
*set
2   # number of sides
3  3
4  3
*elements
*set
4  # number of elements
8  # number of element nodes
1  1  3  11  9  2  7  10  6
2  3  5  13  11  4  8  12  7
3  9  11  19  17  10  15  18  14
4  11  13  21  19  12  16  20  15
# end elements
*nodes
21  # number of nodes
2   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00
    2   0.5000000e+00   0.0000000e+00
    3   1.0000000e+00   0.0000000e+00
    4   1.5000000e+00   0.0000000e+00
    5   2.0000000e+00   0.0000000e+00

    6   0.0000000e+00   0.5000000e+00
    7   1.0000000e+00   0.5000000e+00
    8   2.0000000e+00   0.5000000e+00

    9   0.0000000e+00   1.0000000e+00
   10   0.5000000e+00   1.0000000e+00
   11   1.0000000e+00   1.0000000e+00
   12   1.5000000e+00   1.0000000e+00
   13   2.0000000e+00   1.0000000e+00

   14   0.0000000e+00   1.5000000e+00
   15   1.0000000e+00   1.5000000e+00
   16   2.0000000e+00   1.5000000e+00

   17   0.0000000e+00   2.0000000e+00
   18   0.5000000e+00   2.0000000e+00
   19   1.0000000e+00   2.0000000e+00
   20   1.5000000e+00   2.0000000e+00
   21   2.0000000e+00   2.0000000e+00
