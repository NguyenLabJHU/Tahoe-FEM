*version
1.0
*title
2 cubes with CSE between
*dimensions
16   # number of nodes
3   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   8
2   1   8
5   # number of node sets
# [ID] [nnd]
1   4 # bottom
2   4 # top
3   8 # LL
4   8 # UL
5   4 # top face
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
8   # number of nodes
1 2 3 4 5 6 7 8
*set
8   # number of nodes
9 10 11 12  13 14 15 16
*set
4   # number of nodes
13 14 15 16
# end node sets
*sidesets
*elements
*set
2   # number of elements
8   # number of element nodes
1  1  2  3 4 5 6 7 8
2  9 10 11 12 13 14 15 16  
*set
1   # number of elements
8   # number of element nodes
1 5  6  7  8  9 10 11 12
# end elements
*nodes
16  # number of nodes
3  # number of spatial dimensions
1  0.0  0.0 0.0
2  1.0  0.0 0.0
3  1.0  1.0 0.0
4  0.0  1.0 0.0
5  0.0  0.0 1.0
6  1.0  0.0 1.0
7  1.0  1.0 1.0
8  0.0  1.0 1.0
9  0.0  0.0 1.0
10 1.0  0.0 1.0
11 1.0  1.0 1.0
12 0.0  1.0 1.0
13 0.0  0.0 2.0
14 1.0  0.0 2.0
15 1.1  1.0 2.0
16 0.0  1.0 2.0