*version
1.0
*title
1 blocks
*dimensions
16  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   9   4
2   # number of node sets
# [ID] [nnd]
1 4
2 16
1   # number of side sets
# [ID] [elem set ID] [ns]
1  1  1
# end dimensions
*nodesets
*set
4   # number of nodes
1  4  7 10
*set
16   # number of nodes
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
# end node sets
*sidesets
*set
1
1 2
*elements
*set
9
4
1 1 2 5 4
2 2 3 6 5
3 4 5 8 7
4 5 6 9 8
5 7 8 11 10
6 8 9 12 11
7 9 15 16 12
8 6 14 15 9
9 3 13 14 6
# end elements
*nodes
16
2
1  2.1 -1.0
2  3.1 -1.0
3  4.1 -1.0

4  2.1 0.0
5  3.1 0.0
6  4.1 0.0

7  2.1 1.0
8  3.1 1.0
9  4.1 1.0

10  2.1 2.0
11  3.1 2.0
12  4.1 2.0

13  5.1 -1.0
14  5.1 0.0
15  5.1 1.0
16  5.1 2.0

