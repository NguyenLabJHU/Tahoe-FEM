*version
1.0
*title
one element contact
*dimensions
8 # number of nodes
3 # number of spatial dimensions
1 # number of element sets
# [ID] [nel] [nen]
1 1 8

6 #number of node sets
# [ID] [nnd]
1 4
2 4
3 4
4 4
5 4
6 4

0 # number of side sets
# end dimensions

*nodesets
*set
4 # number of nodes in set 1
1 2 3 4
*set
4 # number of nodes
5 6 7 8
*set
4 # number of nodes 
4 7 6 1
*set
4 # number of nodes in set 4
3 8 7 4
*set
4 # number of nodes in set 5
2 5 8 3
*set
4 # number of nodes
1 6 5 2
# end node sets

*sidesets
*elements
*set
1 # number of elements
8 # number of nodes per element
1 
4 7 8 3 1 6 5 2
# end elements

*nodes
8 #number of nodes
3 # number of spatial dimensions
1  0.5 -0.5 0.5
2  0.5  0.5 0.5
3 -0.5  0.5 0.5
4 -0.5 -0.5 0.5
5  0.5  0.5 -0.5
6  0.5 -0.5 -0.5
7 -0.5 -0.5 -0.5
8 -0.5  0.5 -0.5
