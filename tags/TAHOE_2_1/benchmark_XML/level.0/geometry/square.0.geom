*version
1.0
*title
1 square element
*dimensions
4   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   4
5   # number of node sets
# [ID] [nnd]
1   1
2   1
3   1
4   1
5   4
4   # number of side sets
# [ID] [element set ID] [ns]
1   1   1
2   1   1
3   1   1
4   1   1
# end dimensions
*nodesets
*set
1   # number of nodes
2
*set
1   # number of nodes
3
*set
1   # number of nodes
4
*set
1   # number of nodes
1
*set
4   # numbr of nodes
1 2 3 4
# end node sets
*sidesets
*set
1   # number of sides
# [element] [face]
1  4
*set
1   # number of sides
# [element] [face]
1  1
*set
1   # number of sides
# [element] [face]
1  2
*set
1   # number of sides
# [element] [face]
1  3 
# end side sets
*elements
*set
1   # number of elements
4   # number of element nodes
1  1  2  3  4
# end elements
*nodes
4  # number of nodes
2  # number of spatial dimensions
 1  0.0  0.0
 2  1.0  0.0
 3  1.0  1.0
 4  0.0  1.0
