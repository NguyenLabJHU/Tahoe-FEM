*version
1.0
*title
rectangular grid
*dimensions
25  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1  16  4
4   # number of node sets
# [ID] [nnd]
1  5
2  5
3  5
4  5
0   # number of side sets
# end dimensions
*nodesets
*set
5
1        6       11       16       21
*set
5
21       22       23       24       25
*set
5
5       10       15       20       25
*set
5
1        2        3        4        5
# end node sets
*sidesets
*elements
*set
grid.elems
*nodes
grid.nodes
