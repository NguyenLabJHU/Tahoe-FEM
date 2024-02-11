*version
1.0
*title
4 linear quads with 4 CSE between
*dimensions
16   # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   4   4
2   4   4
11  # number of node sets
# [ID] [nnd]
1   4 # bottom
2   4 # top
3   2 # LU
4   1 # LL
5   2 # RU
6   2 # RL
7   1 # 
8   2 # top-right
9   2 # top-left
10  2 # bottom-right
11  2 # bottom-left
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1 2 3 4
*set
4   # number of nodes
13 14 15 16
*set
2   # number of nodes
9 13
*set
1   # number of nodes
5
*set
2   # number of nodes
12 16 
*set
2   # number of nodes
4 8
*set
1   # number of nodes
1
*set
2   # number of nodes
15 16
*set
2   # number of nodes
13 14
*set
2   # number of nodes
1 2
*set
2   # number of nodes
3 4
# end node sets
*sidesets
*elements
*set
4   # number of elements
4   # number of element nodes
1  1  2   6   5
2  3  4   8   7
3  11 12  16  15
4  9  10  14  13
*set
4   # number of elements
4   # number of element nodes
1  6  2  3   7
2  7  8  12  11
3  14 10 11  15
4  5  6  10  9
# end elements
*nodes
16  # number of nodes
2  # number of spatial dimensions
1   0.0   0.0
2   0.04  0.0
3   0.04  0.0
4   0.1   0.0
5   0.0   0.04
6   0.055 0.07
7   0.055 0.07
8   0.1   0.03
9   0.0   0.04
10  0.055 0.07
11  0.055 0.07
12  0.1   0.03
13  0.0   0.1
14  0.045 0.1
15  0.045 0.1
16  0.1   0.1

