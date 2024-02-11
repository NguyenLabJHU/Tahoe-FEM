*version
1.0
*title
2x2 elements
*dimensions
9  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   4
2   # number of node sets
# [ID] [nnd]
1   3
2   3
0   # number of side sets
# end dimensions

*nodesets
*set
3   # number of nodes
1  2  3  
*set
3   # number of nodes
7  8  9  
# end node sets
*sidesets
*elements
*set
4  # number of elements
4   # number of element nodes
    1          2       5       6       3
    2          1       4       5       2
    3          5       8       9       6
    4          4       7       8       5

# end elements
*nodes
9 # number of nodes
2  # number of spatial dimensions
    1              0              0           
    2              0            0.5              
    3              0              1              
    4              1              0              
    5              1            0.5              
    6              1              1              
    7              2              0              
    8              2            0.5              
    9              2              1              

