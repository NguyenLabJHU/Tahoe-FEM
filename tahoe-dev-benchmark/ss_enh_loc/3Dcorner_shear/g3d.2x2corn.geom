*version
1.0
*title
Geometry File 3D, 2x2x2mesh for corner shear
*dimensions
27 # number of nodes
3 # number of spatial dimensions

1 # number of element sets
# [ID] [nel] [nen]
1 8 8

9 #number of node sets
# [ID] [nnd]
1 7
2 1
3 9
4 9
5 9
6 9
7 9
8 4
9 9

0 # number of side sets
# end dimensions

*nodesets
*set
7 # number of nodes for set 1
1 3 7 9 19 21 25 
*set
1 # number of nodes for set 2
27
*set
9 # number of nodes for set 3
27 24 15 18 3 6 9 21 12
*set
9 # number of nodes for set 4
27 24 23 26 21 20 19 22 25
*set
9 # number of nodes for set 5
1 2 3 4 5 6 7 8 9
*set
9 # number of nodes for set 6
1 4 7 10 13 16 19 22 25
*set
9 # number of nodes for set 7
1 2 3 10 11 12 19 20 21
*set
4 # number of nodes for set 8
17 18 26 27
*set
9 # number of nodes for set 9
17 18 26 27 7 8 9 16 25
# end node sets

*sidesets
*elements
*set
8 # number of elements
8 # number of nodes per element
8 1 2 5 4 10 11 14 13
7 2 3 6 5 11 12 15 14
6 4 5 8 7 13 14 17 16
5 5 6 9 8 14 15 18 17
4 10 11 14 13 19 20 23 22
3 11 12 15 14 20 21 24 23
2 13 14 17 16 22 23 26 25
1 14 15 18 17 23 24 27 26
# end elements
*nodes
27 #number of nodes
3 # number of spatial dimensions
1 0 0 0
2 0.005 0 0
3 0.01 0 0
4 0 0.005 0
5 0.005 0.005 0
6 0.01 0.005 0
7 0 0.01 0
8 0.005 0.01 0
9 0.01 0.01 0
10 0 0 0.005
11 0.005 0 0.005
12 0.01 0 0.005
13 0 0.005 0.005
14 0.005 0.005 0.005
15 0.01 0.005 0.005
16 0 0.01 0.005
17 0.005 0.01 0.005
18 0.01 0.01 0.005
19 0 0 0.01
20 0.005 0 0.01
21 0.01 0 0.01
22 0 0.005 0.01
23 0.005 0.005 0.01
24 0.01 0.005 0.01
25 0 0.01 0.01
26 0.005 0.01 0.01
27 0.01 0.01 0.01



