*version
1.0
*title
2-2x1 blocks 
*dimensions
19  # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen]
1   2   4
2   2   8
2   # number of node sets
# [ID] [nnd]
1   3
2   5
2   # number of side sets
# [ID] [element set ID] [ns]
1   1   2
2   2   2
# end dimensions
*nodesets
*set
3   # number of nodes
1  2  3
*set
5   # number of nodes
15 16 17 18 19
# end node sets
*sidesets
*set
2   # number of sides
# [element] [face]
1  3
2  3
*set
2   # number of sides
# [element] [face]
1  1
2  1
# end side sets
*elements
*set
2   # number of elements
4   # number of element nodes
1  1  2  5  4
2  2  3  6  5
*set
2   # number of elements
8   # number of element nodes
1  7   9  17  15   8  13  16  12
2  9  11  19  17  10  14  18  13
# end elements
*nodes
19 # number of nodes
 2 # number of spatial dimensions
 1  0.0  0.0
 2  1.0  0.0
 3  2.0  0.0

 4  0.0  1.0
 5  1.0  1.0
 6  2.0  1.0

 7  0.0  2.0
 8  0.5  2.0
 9  1.0  2.0
10  1.5  2.0
11  2.0  2.0

12  0.0  2.5
13  1.0  2.5
14  2.0  2.5

15  0.0  3.0
16  0.5  3.0
17  1.0  3.0
18  1.5  3.0
19  2.0  3.0
