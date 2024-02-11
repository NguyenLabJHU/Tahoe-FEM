*version
1.0
*title
Geometry File 3D
*dimensions
12 # number of nodes
3 # number of spatial dimensions
1 # number of element sets
# [ID] [nel] [nen]
1 2 8

4 #number of node sets
# [ID] [nnd]
1 4
2 4
3 2
4 12

0 # number of side sets
# end dimensions

*nodesets
*set
4 # number of nodes in set 1
3 6 9 12
*set
4 # number of nodes
1 4 10 7
*set
2 # number of nodes 
1 7
*set
12 # number of nodes in set 4
1 2 3 4 5 6 7 8 9 10 11 12
# end node sets

*sidesets
*elements
*set
2 # number of elements
8 # number of nodes per element
1 
1 4 5 2 7 10 11 8
2 
2 5 6 3 8 11 12 9
# end elements

*nodes
12 #number of nodes
3 # number of spatial dimensions
1 0 0 0
2 0 0.04 0
3 0 0.08 0
4 0.04 0 0
5 0.04 0.04 0
6 0.04 0.08 0
7 0 0 0.08
8 0 0.04 0.08
9 0 0.08 0.08
10 0.04 0 0.08
11 0.04 0.04 0.08
12 0.04 0.08 0.08
