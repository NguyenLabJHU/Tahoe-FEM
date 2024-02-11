*version
1.0
*title
one q8q4 elements
*dimensions
8   # number of nodes
2   # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   1   8 

5   # number of node sets
# [ID] [nnd]
1   3
2   1
3   3
4   1
5   5



3   # number of side sets
# [ID] [element set ID] [ns]
1   1   1
2   1   1
3   1   1


*nodesets

*set
3  # number of nodes in set 1 (left side nodes)
1 
4
6

*set
1  # number of nodes in set 2 
8

*set
3  # number of nodes in set 3 (right side nodes)
3
5
8

*set
1  # number of nodes in set 4
1

*set
5  # number of nodes in set 5 (left side and bottom nodes)
1
2 
3 
4
6






*sidesets
*set
1   # number of sides
# [element] [face]
1  1
*set
1   # number of sides
# [element] [face]
1  3
*set
1   # number of sides
# [element] [face]
1  4
# end side sets



*elements
*set
1  # number of elements
8  # number of element nodes
1 1 3 8 6 2 5 7 4


*nodes
8  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  0.5e-3  0.0
3  1.0e-3  0.0
4  0.0  0.5e-3
5  1.0e-3  0.5e-3
6  0.0  1.0e-3
7  0.5e-3  1.0e-3
8  1.0e-3  1.0e-3






