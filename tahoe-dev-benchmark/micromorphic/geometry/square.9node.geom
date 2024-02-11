*version

1.0

*title
1 element square with 4 nodes

*dimensions
9  # number of nodes
2   # number of spatial dimensions

1   # number of element sets
# [ID] [nel] [nen]
1   1   9

4   # number of node sets
# [ID] [nnd]
1   3
2   3
3   3
4   4




0   # number of side sets
# end dimensions


*nodesets
*set
3   # number of nodes
1 5 2
*set
3   # number of nodes
4 8 1
*set
3   # number of nodes
3 7 4 
*set
4   # number of nodes
1 2 3 4 
# end node sets



*sidesets


*elements
*set
1  # number of elements
9 # number of element nodes
1 1 2 3 4 5 6 7 8 9 
# end elements


*nodes
9  # number of nodes
2   # number of spatial dimensions
    1   -1.0000000e+00   -1.0000000e+00   
    2    1.0000000e+00   -1.0000000e+00   
    3    1.0000000e+00    1.0000000e+00   
    4   -1.0000000e+00    1.0000000e+00   
    5    0.0000000e+00   -1.0000000e+00  
    6    1.0000000e+00    0.0000000e+00
    7    0.0000000e+00    1.0000000e+00    
    8   -1.0000000e+00    0.0000000e+00
    9    0.0000000e+00    0.0000000e+00

