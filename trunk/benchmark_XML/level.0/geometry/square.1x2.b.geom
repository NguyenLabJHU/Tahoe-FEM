*version
1.0
*title
2 solid elements
*dimensions
8
2
2 # number of element sets
      1      1       4
      2      1       4
3 # number of node sets
     100      2
     101      2
     102      4
1 # number of side sets
# [ID] [element set ID] [ns]
    1    1    4
*nodesets
*set
2
1  2
*set
2
7  8
*set
4
5  6  7  8
*sidesets
*set
4   # number of sides
# [element] [face]
1    1
1    2
1    3
1    4
*elements
*set
1  4
1  1  2  4  3
*set
1  4
1  5  6  8  7
*nodes
8  2
1 -1.0  0.0
2  2.0  0.0
3 -1.0  1.0
4  2.0  1.0

5  0.0  1.5
6  1.0  1.5
7  0.0  2.5
8  1.0  2.5
