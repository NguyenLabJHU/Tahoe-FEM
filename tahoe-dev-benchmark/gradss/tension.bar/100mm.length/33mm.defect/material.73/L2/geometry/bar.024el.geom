*version
1.0
*title
1D bar/rod
*dimensions
25  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   16  2
2   8   2
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
16  # number of elements
2   # number of element nodes
1	2	3
2	3	4
3	4	5
4	5	6
5	6	7
6	7	8
7	8	9
8	9	10
9	18	19
10	19	20
11	20	21
12	21	22
13	22	23
14	23	24
15	24	25
16	25	1
*set
8   # number of elements
2   # number of element nodes
1	10	11
2	11	12
3	12	13
4	13	14
5	14	15
6	15	16
7	16	17
8	17	18
# end elements
*nodes
25  # number of nodes
1   # number of spatial dimensions
2	0
3	4.17
4	8.33
5	12.5
6	16.67
7	20.83
8	25
9	29.17
10	33.33
11	37.5
12	41.67
13	45.83
14	50
15	54.17
16	58.33
17	62.5
18	66.67
19	70.83
20	75
21	79.17
22	83.33
23	87.5
24	91.67
25	95.83
1	100
