*version
1.0
*title
5x5 square grid
*dimensions
36  # number of nodes
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1   25   4
8   # number of node sets
# [ID] [nnd]
1   6
2   6
3   6
4   6
5   1
6   1
7   1
8   1
0   # number of side sets
# end dimensions

*nodesets
*set
6   # number of nodes
1  7  13  19  25  31
*set
6   # number of nodes
31  32  33  34  35  36
*set
6   # number of nodes
6  12  18  24  30  36
*set
6   # number of nodes
1  2  3  4  5  6
*set
1   # number of nodes
1
*set
1   # number of nodes
31
*set
1   # number of nodes
36
*set
1   # number of nodes
6
# end node sets
*sidesets
*elements
*set
25  # number of elements
4   # number of element nodes
       1       1       7       8       2
       2       7      13      14       8
       3      13      19      20      14
       4      19      25      26      20
       5      25      31      32      26
       6       2       8       9       3
       7       8      14      15       9
       8      14      20      21      15
       9      20      26      27      21
      10      26      32      33      27
      11       3       9      10       4
      12       9      15      16      10
      13      15      21      22      16
      14      21      27      28      22
      15      27      33      34      28
      16       4      10      11       5
      17      10      16      17      11
      18      16      22      23      17
      19      22      28      29      23
      20      28      34      35      29
      21       5      11      12       6
      22      11      17      18      12
      23      17      23      24      18
      24      23      29      30      24
      25      29      35      36      30
# end elements
*nodes
36 # number of nodes
2  # number of spatial dimensions
1    -2.5   -1.25
2    -2.5   -0.75
3    -2.5   -0.25
4    -2.5   0.25
5    -2.5   0.75
6    -2.5   1.25
7    -1.5   -1.25
8    -1.5   -0.75
9    -1.5   -0.25
10   -1.5   0.25
11   -1.5   0.75
12   -1.5   1.25
13   -0.5   -1.25
14   -0.5   -0.75
15   -0.5   -0.25
16   -0.5   0.25
17   -0.5   0.75
18   -0.5   1.25
19   0.5    -1.25
20   0.5    -0.75
21   0.5    -0.25
22   0.5    0.25
23   0.5    0.75
24   0.5    1.25
25   1.5    -1.25
26   1.5    -0.75
27   1.5    -0.25
28   1.5    0.25
29   1.5    0.75
30   1.5    1.25
31   2.5    -1.25
32   2.5    -0.75
33   2.5    -0.25
34   2.5    0.25
35   2.5    0.75
36   2.5    1.25