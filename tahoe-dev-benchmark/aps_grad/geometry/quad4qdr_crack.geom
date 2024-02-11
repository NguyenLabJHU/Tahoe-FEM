*version
1.0
*title
2 x 2 patch of quad4 q8q4 elements
*dimensions
23   # number of nodes
2   # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   4   8 

4   # number of node sets
# [ID] [nnd]
1   5
2   10 
3   3
4   3



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
16
19

*set
10  # number of nodes in set 2 (top and bottom nodes) 
1
2
3
4
5
19 
20
21
22
23
 

*set
3  # number of nodes in set 3
5 8 13

*set
3  # number of nodes in set 4
15 18 23




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
8  # number of element nodes
1 1 3 11 9 2 7 10 6
2 3 5 13 11 4 8 12 7 
3 9 11 21 19 10 17 20 16
4 11 15 23 21 14 18 22 17

*nodes
23  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  0.25e-3  0.0
3  0.5e-3  0.0
4  0.75e-3  0.0
5  1.0e-3  0.0
6  0.0  0.25e-3
7  0.5e-3  0.25e-3
8  1.0e-3  0.25e-3
9  0.0  0.5e-3
10  0.25e-3  0.5e-3
11  0.5e-3  0.5e-3
12  0.75e-3  0.5e-3
13  1.0e-3  0.5e-3
14  0.75e-3  0.5e-3
15  1.0e-3  0.5e-3
16  0.0  0.75e-3
17  0.5e-3  0.75e-3
18  1.0e-3  0.75e-3
19  0.0  1.0e-3
20  0.25e-3  1.0e-3
21  0.5e-3  1.0e-3
22  0.75e-3  1.0e-3
23  1.0e-3  1.0e-3





