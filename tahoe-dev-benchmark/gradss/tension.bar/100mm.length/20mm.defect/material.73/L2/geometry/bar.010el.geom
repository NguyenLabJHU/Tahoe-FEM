*version
1.0
*title
1D bar/rod
*dimensions
11  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   8   2
2   2   2
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
8   # number of elements
2   # number of element nodes
1	2	3
2	3	4
3	4	5
4	5	6
5	8	9
6	9	10
7	10	11
8	11	1
*set
2   # number of elements
2   # number of element nodes
1	6	7
2	7	8
# end elements
*nodes
11 # number of nodes
1  # number of spatial dimensions
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
