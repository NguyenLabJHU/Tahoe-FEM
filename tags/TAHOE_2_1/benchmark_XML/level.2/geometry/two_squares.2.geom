*version
1.0
*title
2 squares
*dimensions
8   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   2   4
5   # number of node sets
# [ID] [nnd]
1   2
2   2
3   2
4   2
5   2
2   # number of side sets
# [ID] [assoc elem block] [ns]
1  1  1
2  1  1
# end dimensions
*nodesets
*set
2   # number of nodes
1  2
*set
2   # number of nodes
3  4
*set
2   # number of nodes
5  6
*set
2   # number of nodes
7  8
*set
2   # number of nodes
1  7
# end node sets
*sidesets
*set
1
1 3
*set
1
2 1
*elements
*set
2  # number of elements
4  # number of element nodes
1  1  2  4  3
2  5  6  8  7
# end elements
*nodes
8   # number of nodes
2   # number of spatial dimensions
1 -0.1  0.0
2  5.1  0.0

3 -0.1  1.0
4  5.1  1.0

5  0.0  1.0
6  1.0  1.0

7  0.0  2.0
8  1.0  2.0
