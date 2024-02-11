*version
1.0
*title
1 x 1 element, 5 x 10 plate, 2d bar in tension, no defect
*dimensions
4   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   4
4   # number of node sets
# [ID] [nnd]
1   2
2   2
3   1
3   2
0   # number of side sets
*nodesets
*set
2   # number of nodes
2
3
*set
2   # number of nodes
4
1
*set
1   # number of nodes
2
*set
2   # number of nodes
2
4
# end node sets
*sidesets
*elements
*set
1   # number of elements
4   # number of element nodes
1	2	3	1	4
*nodes
4  # number of nodes
2  # number of spatial dimensions
2	0	0
3	5	0
4	0	10
1	5	10
