*version
1.0
*title
5 x 10 element, 5 x 10 plate, 2d bar in tension, 1 x 1 defect
*dimensions
66  # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   49  4
2   1   4
4   # number of node sets
# [ID] [nnd]
1   6
2   6
3   1
4   11
0   # number of side sets
*nodesets
*set
6   # number of nodes
2
3
4
5
6
7
*set
6   # number of nodes
62
63
64
65
66
1
*set
1   # number of nodes
2
*set
11  # number of nodes
2
8
14
20
26
32
38
44
50
56
62
# end node sets
*sidesets
*elements
*set
49  # number of elements
4   # number of element nodes
1	3	4	10	9
2	4	5	11	10
3	5	6	12	11
4	6	7	13	12
5	8	9	15	14
6	9	10	16	15
7	10	11	17	16
8	11	12	18	17
9	12	13	19	18
10	14	15	21	20
11	15	16	22	21
12	16	17	23	22
13	17	18	24	23
14	18	19	25	24
15	20	21	27	26
16	21	22	28	27
17	22	23	29	28
18	23	24	30	29
19	24	25	31	30
20	26	27	33	32
21	27	28	34	33
22	28	29	35	34
23	29	30	36	35
24	30	31	37	36
25	32	33	39	38
26	33	34	40	39
27	34	35	41	40
28	35	36	42	41
29	36	37	43	42
30	38	39	45	44
31	39	40	46	45
32	40	41	47	46
33	41	42	48	47
34	42	43	49	48
35	44	45	51	50
36	45	46	52	51
37	46	47	53	52
38	47	48	54	53
39	48	49	55	54
40	50	51	57	56
41	51	52	58	57
42	52	53	59	58
43	53	54	60	59
44	54	55	61	60
45	56	57	63	62
46	57	58	64	63
47	58	59	65	64
48	59	60	66	65
49	60	61	1	66
*set
1   # number of elements
4   # number of element nodes
# end elements
1	2	3	9	8
*nodes
66 # number of nodes
2  # number of spatial dimensions
2	0	0
3	1	0
4	2	0
5	3	0
6	4	0
7	5	0
8	0	1
9	1	1
10	2	1
11	3	1
12	4	1
13	5	1
14	0	2
15	1	2
16	2	2
17	3	2
18	4	2
19	5	2
20	0	3
21	1	3
22	2	3
23	3	3
24	4	3
25	5	3
26	0	4
27	1	4
28	2	4
29	3	4
30	4	4
31	5	4
32	0	5
33	1	5
34	2	5
35	3	5
36	4	5
37	5	5
38	0	6
39	1	6
40	2	6
41	3	6
42	4	6
43	5	6
44	0	7
45	1	7
46	2	7
47	3	7
48	4	7
49	5	7
50	0	8
51	1	8
52	2	8
53	3	8
54	4	8
55	5	8
56	0	9
57	1	9
58	2	9
59	3	9
60	4	9
61	5	9
62	0	10
63	1	10
64	2	10
65	3	10
66	4	10
1	5	10
