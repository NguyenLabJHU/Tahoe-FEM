*version
1.0
*title
2 quads with mode-1 CSE (symmetry enforcement replaces element on top. side set gives CSE) 
*dimensions
4   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   4
2   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   1 # LL
1   # number of side sets
# [ID] [elem set ID] [ns]
1 1 1 
# end dimensions
*nodesets
*set
2   # number of nodes
1  2 #3 4
*set
1   # number of nodes
1
# end node sets
*sidesets
*set
1 # 1 side
1 3 # element block 1, side 3
*elements
*set
1   # number of elements
4   # number of element nodes
1  1  2  4  3
# end elements
*nodes
4  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  1.0  0.0
3  0.0  1.0
4  1.0  1.0
