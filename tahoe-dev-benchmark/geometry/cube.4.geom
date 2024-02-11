*version
1.0
*title
1 element cube with 20 nodes
*dimensions
20  # number of nodes
3   # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   1   20

3   # number of node sets
# [ID] [nnd]
1   8
2   8
3   8

0   # number of side sets
# end dimensions

*nodesets
*set
8   # number of nodes
1 2 3 4 9 10 11 12
*set
8   # number of nodes
1 4 8 5 12 20 16 17
*set
8   # number of nodes
1 2 6 5 9 18 13 17
# end node sets

*sidesets

*elements
*set
1  # number of elements
20 # number of element nodes
1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
# end elements

*nodes
20  # number of nodes
3   # number of spatial dimensions
    1   0.0000000e+00   0.0000000e+00   0.0000000e+00
    2   1.0000000e+00   0.0000000e+00   0.0000000e+00
    3   1.0000000e+00   1.0000000e+00   0.0000000e+00
    4   0.0000000e+00   1.0000000e+00   0.0000000e+00
    5   0.0000000e+00   0.0000000e+00   1.0000000e+00
    6   1.0000000e+00   0.0000000e+00   1.0000000e+00
    7   1.0000000e+00   1.0000000e+00   1.0000000e+00
    8   0.0000000e+00   1.0000000e+00   1.0000000e+00
    9   0.5000000e+00   0.0000000e+00   0.0000000e+00
   10   1.0000000e+00   0.5000000e+00   0.0000000e+00
   11   0.5000000e+00   1.0000000e+00   0.0000000e+00
   12   0.0000000e+00   0.5000000e+00   0.0000000e+00
   13   0.5000000e+00   0.0000000e+00   1.0000000e+00
   14   1.0000000e+00   0.5000000e+00   1.0000000e+00
   15   0.5000000e+00   1.0000000e+00   1.0000000e+00
   16   0.0000000e+00   0.5000000e+00   1.0000000e+00
   17   0.0000000e+00   0.0000000e+00   0.5000000e+00
   18   1.0000000e+00   0.0000000e+00   0.5000000e+00
   19   1.0000000e+00   1.0000000e+00   0.5000000e+00
   20   0.0000000e+00   1.0000000e+00   0.5000000e+00
