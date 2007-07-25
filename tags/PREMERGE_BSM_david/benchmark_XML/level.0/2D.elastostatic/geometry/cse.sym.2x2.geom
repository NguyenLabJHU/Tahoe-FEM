*version
1.0
*title
2x2 quads with mode-1 CSE (symmetry enforcement replaces element on top. side set gives CSE) 
*dimensions
9   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   4
2   # number of node sets
# [ID] [nnd]
1   3 # bottom
2   1 # LL
1   # number of side sets
# [ID] [elem set ID] [ns]
1 1 2
# end dimensions
*nodesets
*set
3   # number of nodes
1  2  3
*set
1   # number of nodes
1
# end node sets
*sidesets
*set
2
3 3
4 3
*elements
*set
4   # number of elements
4   # number of element nodes
1  1  2  5  4
2  2  3  6  5
3  4  5  8  7
4  5  6  9  8
# end elements
*nodes
9  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  2.0  0.0

4  0.0  1.0
5  1.0  1.0
6  2.0  1.0

7  0.0  2.0
8  1.0  2.0
9  2.0  2.0
