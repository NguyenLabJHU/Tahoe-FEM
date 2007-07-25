*version
1.0
*title
2 blocks
*dimensions
32  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   6   8
4   # number of node sets
# [ID] [nnd]
1   12
2   12
3   4
4   4
2   # number of side sets
# [ID] [assoc elem block] [ns]
1  1  5
2  1  1
# end dimensions
*nodesets
*set
12   # number of nodes
1  2  3  4 5  6  7  8 9 10 11 12
*set
12   # number of nodes
13 14 15 16 17 18 19 20 21 22 23 24
*set
4   # number of nodes
25 26 27 28
*set
4   # number of nodes
29 30 31 32
# end node sets
*sidesets
*set
5
1 2
2 2
3 2
4 2
5 2
*set
1
6 1
*elements
*set
6  # number of elements
8  # number of element nodes
1  1  2  4  3  13  14  16  15
2  3  4  6  5  15  16  18  17
3  5  6  8  7  17  18  20  19
4  7  8  10 9  19  20  22  21
5  9  10  12  11  21  22  24 23  
6  25 26 27 28 29 30 31 32
# end elements
*nodes
32  # number of nodes
3   # number of spatial dimensions
1  -0.1 -0.1  0.0
2   1.1 -0.1  0.0
3  -0.1  0.94 0.0
4   1.1  0.94 0.0
5  -0.1  1.98 0.0
6   1.1  1.98 0.0
7  -0.1  3.02 0.0
8   1.1  3.02 0.0
9  -0.1  4.06 0.0
10  1.1  4.06 0.0
11 -0.1  5.1  0.0
12  1.1  5.1  0.0

13 -0.1 -0.1  1.0
14  1.1 -0.1  1.0
15 -0.1  0.94 1.0
16  1.1  0.94 1.0
17 -0.1  1.98 1.0
18  1.1  1.98 1.0
19 -0.1  3.02 1.0
20  1.1  3.02 1.0
21 -0.1  4.06 1.0
22  1.1  4.06 1.0
23 -0.1  5.1  1.0
24  1.1  5.1  1.0


25  0.0  0.0  1.0
26  1.0  0.0  1.0
27  1.0  1.0  1.0
28  0.0  1.0  1.0

29  0.0  0.0  2.0
30  1.0  0.0  2.0
31  1.0  1.0  2.0
32  0.0  1.0  2.0
