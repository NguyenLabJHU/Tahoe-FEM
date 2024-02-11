*version
1.0
*title
2 x 2 x 2 element cube patch
*dimensions
27  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   8   8

13   # number of node sets
# [ID] [nnd]
1   9
2   9
3   9
4   9
5   27
6   18
7   9
8   9
9   26
10  19
11  15
12  12
13  19


0   # number of side sets
# end dimensions
*nodesets
*set
9   # number of nodes in set 1
1  2  3  
4  5  6
7  8  9
*set
9   # number of nodes in set 2
19 20 21
22 23 24
25 26 27
*set
9   # number of nodes in set 3
 1  2  3
10 11 12
19 20 21
*set
9   # number of nodes in set 4
 1  4  7
10 13 16
19 22 25
*set
27  # number of nodes in set 5
1  2  3  
4  5  6
7  8  9
10 11 12
13 14 15
16 17 18
19 20 21
22 23 24
25 26 27
*set
18  # number of nodes in set 6
10 11 12
13 14 15
16 17 18
19 20 21
22 23 24
25 26 27
*set
9  # number of nodes in set 7
 7  8  9
16 17 18
25 26 27
*set
9  # number of nodes in set 8
 3  6  9
12 15 18
21 24 27
*set
26  # number of nodes in set 9
1  2  3  
4  5  6
7  8  9
10 11 12
13 15
16 17 18
19 20 21
22 23 24
25 26 27
*set
19  # number of nodes in set 10
1  2  3  
4  5  6
7  8  9
10 11 12
19 20 21
13 16
22 25
*set
15  # number of nodes in set 11
21
24
25 26 27
 7  8  9
16 17 18
 3  6
12 15
*set
12  # number of nodes in set 12
 1  2
 4  5
10 11
13 14
19 20
22 23
*set
19  # number of nodes in set 13
3 6 9 8 7
12 15 18 17 16
21 24 27 26 25
19 20 22 23
# end node sets

*sidesets
*elements
*set
cube.1.elem
# end elements
*nodes
cube.1.node
