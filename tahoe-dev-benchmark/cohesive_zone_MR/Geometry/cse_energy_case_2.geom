*version
1.0
*title
5 quad elements with 8 CSEs between
*dimensions
20  # number of nodes
2   # number of spatial dimensions
4   # number of element sets
# [ID] [nel] [nen]
1   1   4
2   4   4
3   4   4
4   4   4
6   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   2 # LL
4   2 # RL
5   1 # LBN
6   1 # RBN
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
12  11
*set
2   # number of nodes
13  16
*set
2   # number of nodes
6  7
*set
1   # number of nodes
1
*set
1   # number of nodes
2
# end node sets
*sidesets
*elements
*set
1   # number of elements
4   # number of element nodes
1  17  18  19  20
*set
4   # number of elements
4   # number of element nodes
1  1  2  3  4
2  5  6  7  8
3  9  10  11  12
4  13  14  15  16
*set
4   # number of elements
4   # number of element nodes
1  1  4  14  13
2  3  2  6  5
3  8  7  11  10
4  16  15  9  12
*set
4   # number of elements
4   # number of element nodes
1  15  14  17  20
2  4  3  18  17
3  19  18  5  8
4  20  19  10  9
# end elements
*nodes
20  # number of nodes
2  # number of spatial dimensions
1  0.00  0.00
2  1.00  0.00
3  0.75  0.30
4  0.25  0.50
5  0.75  0.30
6  1.00  0.00
7  1.00  1.00
8  0.70  0.75
9  0.20  0.80
10  0.70  0.75
11  1.00  1.00
12  0.00  1.00
13  0.00  0.00
14  0.25  0.50
15  0.20  0.80
16  0.00  1.00
17  0.25  0.50
18  0.75  0.30
19  0.70  0.75
20  0.20  0.80
