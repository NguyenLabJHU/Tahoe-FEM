*version
1.0
*title
1D bar/rod
*dimensions
13  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   8   2
2   4   2
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
5	10	11
6	11	12
7	12	13
8	13	1
*set
4   # number of elements
2   # number of element nodes
1	6	7
2	7	8
3	8	9
4	9	10
# end elements
*nodes
13  # number of nodes
1   # number of spatial dimensions
2	0
3	8.33
4	16.67
5	25
6	33.33
7	41.67
8	50
9	58.33
10	66.67
11	75
12	83.33
13	91.67
1	100
