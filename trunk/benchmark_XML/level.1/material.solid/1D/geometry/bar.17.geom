*version
1.0
*title
1D bar/rod
*dimensions
18  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   16  2
2   1   2
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
9	11	12
10	12	13
11	13	14
12	14	15
13	15	16
14	16	17
15	17	18
16	18	1
*set
1   # number of elements
2   # number of element nodes
1  10  11
# end elements
*nodes
18  # number of nodes
1   # number of spatial dimensions
2	0.000000
3	0.294118
4	0.588235
5	0.882353
6	1.176471
7	1.470588
8	1.764706
9	2.058824
10	2.352941
11	2.647059
12	2.941176
13	3.235294
14	3.529412
15	3.823529
16	4.117647
17	4.411765
18	4.705882
1	5.000000
