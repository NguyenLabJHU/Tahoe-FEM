*version
1.0
*title
4 x 4 square in 3 blocks
*dimensions
33  # number of nodes
2   # number of spatial dimensions
4   # number of element sets
# [ID] [nel] [nen]
1  8  4
2  4  4
3  4  4
4  6  4
0   # number of node sets
# [ID] [nnd]
1   # number of side sets
# [ID] [element set ID] [ns]
1   3   2
# end dimensions
*nodesets
# end node sets
*sidesets
*set
2   # number of sides
# [element] [face]
3  3
4  3
# end side sets
*elements
*set
8   # number of elements
4   # number of element nodes
1  1  2  7  6
2  2  3  8  7
3  3  4  9  8
4  4  5  10 9

5  6  7  12  11
6  7  8  13  12
7  8  9  14  13
8  9  10  15   14
*set
4   # number of elements
4   # number of element nodes
1  16  17  23  22
2  17  18  24  23
3  22  23  29  28
4  23  24  30  29
*set
4   # number of elements
4   # number of element nodes
1  19  20  26  25
2  20  21  27  26
3  25  26  32  31
4  26  27  33  32
*set
6   # number of elements
4   # number of element nodes
1  11  12  17  16
2  12  13  18  17
3  13  14  20  19
4  14  15  21  20
5  24  18  19  25
6  30  24  25  31
# end elements
*nodes
33 # number of nodes
2  # number of spatial dimensions
 1  0.0  0.0
 2  1.0  0.0
 3  2.0  0.0
 4  3.0  0.0
 5  4.0  0.0

 6  0.0  1.0
 7  1.0  1.0
 8  2.0  1.0
 9  3.0  1.0
10  4.0  1.0

11  0.0  2.0
12  1.0  2.0
13  2.0  2.0
14  3.0  2.0
15  4.0  2.0

16  0.0  2.0
17  1.0  2.0
18  2.0  2.0
19  2.0  2.0
20  3.0  2.0
21  4.0  2.0

22  0.0  3.0
23  1.0  3.0
24  2.0  3.0
25  2.0  3.0
26  3.0  3.0
27  4.0  3.0

28  0.0  4.0
29  1.0  4.0
30  2.0  4.0
31  2.0  4.0
32  3.0  4.0
33  4.0  4.0
