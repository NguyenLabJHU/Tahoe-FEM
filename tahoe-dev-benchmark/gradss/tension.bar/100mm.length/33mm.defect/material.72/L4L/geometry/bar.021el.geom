*version
1.0
*title
1D bar/rod
*dimensions
22  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   14  2
2   7   2
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
14  # number of elements
2   # number of element nodes
1	2	3
2	3	4
3	4	5
4	5	6
5	6	7
6	7	8
7	8	9
8	16	17
9	17	18
10	18	19
11	19	20
12	20	21
13	21	22
14	22	1
*set
7   # number of elements
2   # number of element nodes
1	9	10
2	10	11
3	11	12
4	12	13
5	13	14
6	14	15
7	15	16
# end elements
*nodes
22  # number of nodes
1   # number of spatial dimensions
2	0.00
3	4.76
4	9.52
5	14.29
6	19.05
7	23.81
8	28.57
9	33.33
10	38.10
11	42.86
12	47.62
13	52.38
14	57.14
15	61.90
16	66.67
17	71.43
18	76.19
19	80.95
20	85.71
21	90.48
22	95.24
1	100.00
