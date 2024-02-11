*version
1.0
*title
1D bar/rod
*dimensions
19  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   8   3
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
8   # number of elements
3   # number of element nodes
1	2	4	3
2	4	6	5
3	6	8	7
4	8	10	9
5	12	14	13
6	14	16	15
7	16	18	17
8	18	1	19
*set
1   # number of elements
3   # number of element nodes
1	10	12	11
# end elements
*nodes
19  # number of nodes
1   # number of spatial dimensions
2	0
3	5.555555556
4	11.11111111
5	16.66666667
6	22.22222222
7	27.77777778
8	33.33333333
9	38.88888889
10	44.44444444
11	50
12	55.55555556
13	61.11111111
14	66.66666667
15	72.22222222
16	77.77777778
17	83.33333333
18	88.88888889
19	94.44444444
1	100
