*version
1.0
*title
2 beams with surface layer between
*dimensions
30  # number of nodes
2   # number of spatial dimensions
3   # number of element sets
# [ID] [nel] [nen]
1   8   4
2   8   4
3   3   4
2   # number of node sets
# [ID] [nnd]
1   9
2   9
0   # number of side sets
# end dimensions
*nodesets
*set
9   # number of nodes
1  2  3  4  5  6  10  11  15
*set
9   # number of nodes
16 20 21 25 26 27 28 29 30
# end node sets
*sidesets
*elements
*set
8  # number of elements
4  # number of element nodes
1  1  2  7  6
2  2  3  8  7
3  3  4  9  8
4  4  5 10  9
5  6  7 12 11
6  7  8 13 12
7  8  9 14 13
8  9 10 15 14
*set
8  # number of elements
4  # number of element nodes
1  16 17 22 21
2  17 18 23 22
3  18 19 24 23
4  19 20 25 24
5  21 22 27 26
6  22 23 28 27
7  23 24 29 28
8  24 25 30 29
*set
3  # number of elements
4  # number of element nodes
1 12 13 18 17
2 13 14 19 18
3 14 15 20 19
# end elements
*nodes
30  # number of nodes
2   # number of spatial dimensions
    1  -2.0000000e+00  -2.0000000e+00
    2  -1.0000000e+00  -2.0000000e+00
    3   0.0000000e+00  -2.0000000e+00
    4   1.0000000e+00  -2.0000000e+00
    5   2.0000000e+00  -2.0000000e+00
    6  -2.0000000e+00  -1.0000000e+00
    7  -1.0000000e+00  -1.0000000e+00
    8   0.0000000e+00  -1.0000000e+00
    9   1.0000000e+00  -1.0000000e+00
   10   2.0000000e+00  -1.0000000e+00
   11  -2.0000000e+00  -1.0000000e-12
   12  -1.0000000e+00  -1.0000000e-12
   13   0.0000000e+00  -1.0000000e-12
   14   1.0000000e+00  -1.0000000e-12
   15   2.0000000e+00  -1.0000000e-12
   16  -2.0000000e+00   1.0000000e-12
   17  -1.0000000e+00   1.0000000e-12
   18   0.0000000e+00   1.0000000e-12
   19   1.0000000e+00   1.0000000e-12
   20   2.0000000e+00   1.0000000e-12
   21  -2.0000000e+00   1.0000000e+00
   22  -1.0000000e+00   1.0000000e+00
   23   0.0000000e+00   1.0000000e+00
   24   1.0000000e+00   1.0000000e+00
   25   2.0000000e+00   1.0000000e+00
   26  -2.0000000e+00   2.0000000e+00
   27  -1.0000000e+00   2.0000000e+00
   28   0.0000000e+00   2.0000000e+00
   29   1.0000000e+00   2.0000000e+00
   30   2.0000000e+00   2.0000000e+00
