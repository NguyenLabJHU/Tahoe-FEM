*version
1.0
*title
1D bar/rod
*dimensions
43  # number of nodes
1   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   14  3
2   7   3
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
3   # number of element nodes
1	2	4	3
2	4	6	5
3	6	8	7
4	8	10	9
5	10	12	11
6	12	14	13
7	14	16	15
8	30	32	31
9	32	34	33
10	34	36	35
11	36	38	37
12	38	40	39
13	40	42	41
14	42	1	43
*set
7   # number of elements
3   # number of element nodes
1	16	18	17
2	18	20	19
3	20	22	21
4	22	24	23
5	24	26	25
6	26	28	27
7	28	30	29
# end elements
*nodes
43  # number of nodes
1   # number of spatial dimensions
2	0.00
3	2.38
4	4.76
5	7.14
6	9.52
7	11.90
8	14.29
9	16.67
10	19.05
11	21.43
12	23.81
13	26.19
14	28.57
15	30.95
16	33.33
17	35.71
18	38.10
19	40.48
20	42.86
21	45.24
22	47.62
23	50.00
24	52.38
25	54.76
26	57.14
27	59.52
28	61.90
29	64.29
30	66.67
31	69.05
32	71.43
33	73.81
34	76.19
35	78.57
36	80.95
37	83.33
38	85.71
39	88.10
40	90.48
41	92.86
42	95.24
43	97.62
1	100.00
