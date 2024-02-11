*version
1.0
*title
1x1 square grid
*dimensions
8  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   8
3   # number of node sets
# [ID] [nnd]
1   2
2   1
3   3
0   # number of side sets
# end dimensions

*nodesets
*set
2   # number of nodes: fixed in y direction
2  4
*set
1   # number of nodes: fixed in x and y direction
1  
*set
3   # number of nodes: displacement applied
6  7  8
# end node sets
*sidesets
*elements
*set
1  # number of elements
8   # number of element nodes
       1       8       6       1       4       7       3       2       5
# end elements
*nodes
8 # number of nodes
2  # number of spatial dimensions
    1               0               0               
    2          0.0375               0               
    3               0           0.075               
    4           0.075               0               
    5           0.075           0.075               
    6               0            0.15               
    7          0.0375            0.15               
    8           0.075            0.15               

