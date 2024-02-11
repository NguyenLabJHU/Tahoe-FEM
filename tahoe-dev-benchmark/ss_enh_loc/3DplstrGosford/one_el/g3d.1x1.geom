*version
1.0
*title
Gosford Sandstone Geometry File 3D
*dimensions
8 # number of nodes
3 # number of spatial dimensions
1 # number of element sets
# [ID] [nel] [nen]
1 1 8

17 #number of node sets
# [ID] [nnd]
1 4
2 4
3 4
4 4
5 8
6 4
7 4
8 2
9 2
10 1
11 1
12 1
13 1
14 1
15 1
16 1
17 2

0 # number of side sets
# end dimensions

*nodesets
*set
4 # number of nodes in set 1
1 3 5 7
*set
4 # number of nodes
2 4 6 8
*set
4 # number of nodes 
1 2 5 6
*set
4 # number of nodes in set 4
3 4 7 8
*set
8 # number of nodes in set 5
1 2 3 4 5 6 7 8 
*set
4 # number of nodes
1 2 3 4
*set
4 # number of nodes
5 6 7 8
*set
2 # number of nodes in set 8
1 5
*set
2 # number of nodes
2 6
*set
1 # number of nodes
1
*set
1 # number of nodes
2
*set
1 # number of nodes
3
*set
1 # number of nodes
4
*set
1 # number of nodes
5
*set
1 # number of nodes
6
*set
1 # number of nodes
7
*set
2 # number of nodes
1 5
# end node sets

*sidesets
*elements
*set
1 # number of elements
8 # number of nodes per element
1 
1 3 4 2 5 7 8 6
# end elements

*nodes
8 #number of nodes
3 # number of spatial dimensions
1 0 0 0
2 0 0.08 0
3 0.04 0 0
4 0.04 0.08 0
5 0 0 0.08
6 0 0.08 0.08
7 0.04 0 0.08
8 0.04 0.08 0.08
