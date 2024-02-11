*version
1.0
*title
simple shear in 3D
*dimensions
16 # number of nodes
3 # number of spatial dimensions
2 # number of element sets
# [ID] [nel] [nen]
1 1 8
2 2 8

3 #number of node sets
# [ID] [nnd]
1 16
2 4
3 4

0 # number of side sets
# end dimensions

*nodesets
*set
16 # number of nodes in set 1
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
*set
4 # number of nodes in set 2
7 8 15 16
*set
4 # number of nodes in set 3
1 2 9 10
# end node sets

*sidesets

*elements
*set
1 # number of elements
8 # number of nodes per element
1  3 4 6 5 11 12 14 13
*set
2 # number of elements
8 # number of nodes per element
1  1 2 4 3  9 10 12 11
2  5 6 8 7 13 14 16 15
# end elements

*nodes
16 # number of nodes
3  # number of spatial dimensions
1  0 0 0
2  1 0 0
3  0 1 0
4  1 1 0
5  0 2 0
6  1 2 0
7  0 3 0
8  1 3 0
9  0 0 1
10 1 0 1
11 0 1 1
12 1 1 1
13 0 2 1
14 1 2 1
15 0 3 1
16 1 3 1



