*version
1.0
*title
2 x 2 patch of q9q4 elements
*dimensions
25   # number of nodes
2   # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   4   9

10   # number of node sets
# [ID] [nnd]
1   5
2   10 
3   5
4   9
5   1
6   5
7   5
8   1
9   1
10  9





3   # number of side sets
# [ID] [element set ID] [ns]
1   1   2
2   1   2
3   1   2


*nodesets

*set
5  # number of nodes in set 1 (left side nodes)
1 
6
9
14
17

*set
10  # number of nodes in set 2 (top and bottom nodes) 
1
2
3
4
5
17 
18 
19 
20
21
 
*set
5  # number of nodes in set 3 (right side nodes)
5
8
13
16
21


*set
9  # number of nodes in set 4 
1 
6
9
14
17
18
19
20
21


*set
1  # number of nodes in set 5 
5


*set
5  # number of nodes in set 6 
17
18
19
20
21


*set
5  # number of nodes in set 7
1
2
3
4
5

*set
1  # number of nodes in set 8
17

*set
1  # number of nodes in set 9
21

*set
9  # number of nodes in set 10
1
2
3
4
5
6
9
14
17






*sidesets
*set
2   # number of sides
# [element] [face]
1  1
2  1
*set
2   # number of sides
# [element] [face]
3  3
4  3
*set
2   # number of sides
# [element] [face]
1  4
3  4
# end side sets



*elements
*set
4  # number of elements
9  # number of element nodes
1 1 3 11 9 2 7 10 6 22
2 3 5 13 11 4 8 12 7 23
3 9 11 19 17 10 15 18 14 24
4 11 13 21 19 12 16 20 15 25

*nodes
25  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  0.2e-3  0.0
3  0.4e-3  0.0
4  0.7e-3  0.0
5  1.0e-3  0.0
6  0.0  0.2e-3
7  0.5e-3  0.3e-3
8  1.0e-3  0.25e-3
9  0.0  0.4e-3
10  0.3e-3  0.5e-3
11  0.6e-3  0.6e-3
12  0.8e-3  0.55e-3
13  1.0e-3  0.5e-3
14  0.0  0.7e-3
15  0.55e-3  0.8e-3
16  1.0e-3  0.75e-3
17  0.0  1.0e-3
18  0.25e-3  1.0e-3
19  0.5e-3  1.0e-3
20  0.75e-3  1.0e-3
21  1.0e-3  1.0e-3
22  0.25e-3 0.25e-3
23  0.75e-3 0.25e-3
24  0.25e-3 0.75e-3
25  0.75e-3 0.75e-3





