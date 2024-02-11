*version
1.0
*title
1D bar/rod
*dimensions
11  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   4   3
2   1   3
2   # number of node sets
# [ID] [nnd]
1   1
2   1
0   # number of side sets
*nodesets
*set
1   # number of nodes
2
*set
1   # number of nodes
1
# end node sets
*sidesets
*elements
*set
4   # number of elements
3   # number of element nodes
1	2	4	3
2	4	6	5
3	8	10	9
4	10	1	11
*set
1   # number of elements
3   # number of element nodes
1	6	8	7
# end elements
*nodes
11  # number of nodes
1   # number of spatial dimensions
2	0
3	10
4	20
5	30
6	40
7	50
8	60
9	70
10	80
11	90
1	100
