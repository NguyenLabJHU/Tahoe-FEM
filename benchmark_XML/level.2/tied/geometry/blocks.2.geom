*version
1.0
*title
3 x 4 rectangular region in 3 x 2 blocks
*dimensions
20 # number of nodes
 2 # number of spatial dimensions
2 # number of element sets
# [ID]  [nel]  [nen]
   1     6      4
   2     6      4
4 # number of node sets
# [ID]  [nd]
   1     4
   2     4
   3     4
   4     4
0 # number of side sets
*nodesets
*set
4
1  2  3  4
*set
4
17 18 19 20
*set
4
1  5  9  13  17
*set
4
4  8  12  16  20
*sidesets
*elements
*set
6  4
1  1  2  6  5
2  2  3  7  6
3  3  4  8  7
4  5  6  10  9
5  6  7  11  10
6  7  8  12  11
*set
6  4
1  9  10  14  13
2  10  11  15  14
3  11  12  16  15
4  13  14  18  17
5  14  15  19  18
6  15  16  20  19
# end node sets
*nodes
20 2
1  0.0  0.0
2  1.0  0.0
3  2.0  0.0
4  3.0  0.0

5  0.0  1.0
6  1.0  1.0
7  2.0  1.0
8  3.0  1.0

 9  0.0  2.0
10  1.0  2.0
11  2.0  2.0
12  3.0  2.0

13  0.0  3.0
14  1.0  3.0
15  2.0  3.0
16  3.0  3.0

17  0.0  4.0
18  1.0  4.0
19  2.0  4.0
20  3.0  4.0
