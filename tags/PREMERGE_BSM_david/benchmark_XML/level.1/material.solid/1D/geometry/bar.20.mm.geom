*version
1.0
*title
1D bar/rod
*dimensions
21  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   16  2
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
9	14	15
10	15	16
11	16	17
12	17	18
13	18	19
14	19	20
15	20	21
16	21	1
*set
4   # number of elements
2   # number of element nodes
1       10      11
2       11      12
3       12      13
4	13	14
# end elements
*nodes
21  # number of nodes
1   # number of spatial dimensions
2	0.000
3	5.000
4	10.000
5	15.000
6	20.000
7	25.000
8	30.000
9	35.000
10	40.000
11	45.000
12	50.000
13	55.000
14	60.000
15	65.000
16	70.000
17	75.000
18	80.000
19	85.000
20	90.000
21	95.000
1	100.000
