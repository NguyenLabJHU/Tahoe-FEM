*version
1.0
*title
2 cubes
*dimensions
16  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   8
4   # number of node sets
# [ID] [nnd]
1   4
2   4
3   4
4   4
2   # number of side sets
# [ID] [assoc elem block] [ns]
1  1  1
2  1  1
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  3  4
*set
4   # number of nodes
5  6  7  8
*set
4   # number of nodes
9 10 11 12
*set
4   # number of nodes
13 14 15 16
# end node sets
*sidesets
*set
1
1 2
*set
1
2 1
*elements
*set
2  # number of elements
8  # number of element nodes
1  1  2  3  4  5  6  7  8
2  9 10 11 12 13 14 15 16
# end elements
*nodes
16  # number of nodes
3   # number of spatial dimensions
1  0.0  0.0  0.0
2  1.0  0.0  0.0
3  1.0  1.0  0.0
4  0.0  1.0  0.0

5  0.0  0.0  1.0
6  1.0  0.0  1.0
7  1.0  1.0  1.0
8  0.0  1.0  1.0

9   0.0  0.0  1.0
10  1.0  0.0  1.0
11  1.0  1.0  1.0
12  0.0  1.0  1.0

13  0.0  0.0  2.0
14  1.0  0.0  2.0
15  1.0  1.0  2.0
16  0.0  1.0  2.0
