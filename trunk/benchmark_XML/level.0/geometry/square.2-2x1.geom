*version
1.0
*title
2-2x1 blocks 
*dimensions
12  # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   4
2   2   4
2   # number of node sets
# [ID] [nnd]
1   3
2   3
2   # number of side sets
# [ID] [element set ID] [ns]
1   1   2
2   2   2
# end dimensions
*nodesets
*set
3   # number of nodes
1  2  3
*set
3   # number of nodes
10 11 12
# end node sets
*sidesets
*set
2   # number of sides
# [element] [face]
1  3
2  3
*set
2   # number of sides
# [element] [face]
1  1
2  1
# end side sets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  5  4
2  2  3  6  5
*set
2   # number of elements
4   # number of element nodes
1  7  8  11  10
2  8  9  12  11
# end elements
*nodes
12 # number of nodes
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

10  0.0  3.0
11  1.0  3.0
12  2.0  3.0
