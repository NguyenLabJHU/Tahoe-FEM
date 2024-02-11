*version
1.0
*title
1 x 1 x 1 element patch
*dimensions
8   # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   1   8
4   # number of node sets
# [ID] [nnd]
1   4
2   4
3   4
4   1
0   # number of side sets
# end dimensions
*nodesets
*set
4   # number of nodes
1  2  4  3  
5  6  8  7
*set
4   # number of nodes
3  4  7  8
*set
4   # number of nodes
1  2  5  6
*set
1   # number of nodes
1
# end node sets
*sidesets
*elements
*set
bar.001elx001elx001el.elem
# end elements
*nodes
bar.001elx001elx001el.node
