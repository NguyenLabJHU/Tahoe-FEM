*version
1.0
*title
one linear quad
*dimensions
4   # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   4

5   # number of node sets
# [ID] [nnd]
1   2 # bottom
2   2 # top
3   2
4   1 
5   1 


0   # number of side sets

*nodesets
*set
2   # number of nodes
1 2
*set
2   # number of nodes
3 4
*set
2   # number of nodes
1 4
*set
1   # number of nodes
2
*set
1   # number of nodes
3

*sidesets
*elements
*set
1   # number of elements
4   # number of element nodes
1 1 2 3 4
*nodes
4  # number of nodes
2  # number of spatial dimensions
1  0.0  0.0
2  0.04  0.0
3  0.04  0.08
4  0.0  0.08
