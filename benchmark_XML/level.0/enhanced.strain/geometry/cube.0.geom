*version
1.0
*title
4 square element patch
*dimensions
18  # number of nodes
3   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   4   8
2   # number of node sets
# [ID] [nnd]
1   16
2   2
0   # number of side sets
# end dimensions
*nodesets
*set
16  # number of nodes
1 2 3 4 
6 7 8 9
10 11 12 13
15 16 17 18
*set
2  # number of nodes
5 14
# end node sets
*sidesets
*elements
*set
4   # number of elements
8   # number of element nodes
1  1  2  5  4  10  11  14  13
2  2  3  6  5  11  12  15  14
3  4  5  8  7  13  14  17  16
4  5  6  9  8  14  15  18  17

# end elements
*nodes
18 # number of nodes
3  # number of spatial dimensions
1  0.0e-02  0.0e-00  0.0e-00
2  1.0e-02  0.0e-00  0.0e-00
3  2.0e-02  0.0e-00  0.0e-00

4  0.0e-02  1.0e-00  0.0e-00
5  1.0e-02  1.0e-00  0.0e-00
6  2.0e-02  1.0e-00  0.0e-00

7  0.0e-02  2.0e-00  0.0e-00
8  1.0e-02  2.0e-00  0.0e-00
9  2.0e-02  2.0e-00  0.0e-00

10 0.0e-02  0.0e-00  1.0e-00
11 1.0e-02  0.0e-00  1.0e-00
12 2.0e-02  0.0e-00  1.0e-00

13 0.0e-02  1.0e-00  1.0e-00
14 1.0e-02  1.0e-00  1.0e-00
15 2.0e-02  1.0e-00  1.0e-00

16 0.0e-02  2.0e-00  1.0e-00
17 1.0e-02  2.0e-00  1.0e-00
18 2.0e-02  2.0e-00  1.0e-00
