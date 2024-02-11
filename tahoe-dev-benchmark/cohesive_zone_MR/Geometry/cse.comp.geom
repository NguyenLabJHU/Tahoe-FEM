*version
1.0
*title
2 linear quads with CSE between
*dimensions
8   # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   4
2   1   4
11   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   1 # LL
4   1 # UL
5   2 # followers
6   2 # leaders
7   2
8   2
9   2
10  1
11  1
0   # number of side sets
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
7  8
*set
1   # number of nodes
1
*set
1   # number of nodes
7
*set
2   # number of nodes
3 5 
*set
2   # number of nodes
5 6
*set
2   # number of nodes
1 3
*set
2   # number of nodes
5 7
*set
2   # number of nodes
6 8 
*set
1   # number of nodes
4
*set
1   # number of nodes
2
# end node sets
*sidesets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  4  3
2  5  6  8  7
*set
1   # number of elements
4   # number of element nodes
1  3  4  6  5
# end elements
*nodes
8  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  0.4  0.0
3  0.0  1.5
4  0.4  0.5
5  0.0  1.5
6  0.4  0.5
7  0.0  2.0
8  0.4  2.0
