*version 1.0
*title
ring between platens
*dimensions
707    # number of nodes
2      # spatial dimensions
2      # number of element blocks
1  320  4
2  235  4
6      # number of node sets
1  1
2  5
3  6
4  10
5  14
6  54
4      # number of side sets
1  1  12
2  1  12
3  2  13
4  2  13
*nodesets
*set
1
410
*set
5
410  369  328  287  246
*set
6
165  462  465  464  463  412
*set
10
1    42   83  124  165
412  462  463 464  465
*set
14
466  529  528  527  526  525  
524  523  522  521  520  519  
518  517
*set
54
370  371  372  373  374  375  
376  377  378  379  380  381  
382  1    2    3    4    5    
6    7    8    9   10   11   
12   13  462  467  468  469  
470  471  472  473  474  475  
476  477  478  479    
466  529  528  527  526  
525  524  523  522  521  
520  519  518  517
*sidesets
*set
12
281  3  282  3
283  3  284  3
285  3  286  3
287  3  288  3
289  3  290  3
291  3  292  3
*set
12
1   1   2  1
3   1   4  1
5   1   6  1
7   1   8  1
9   1  10  1
11  1  12  1
*set
13
64  1  63  4
62  4  61  4
60  4  59  4
58  4  57  4
56  4  55  4
54  4  53  4
52  4
*set
13
1   4   2  4
3   4   4  4
5   4   6  4
7   4   8  4
9   4  10  4
1   4  12  4
13  4
*nodes
ring.geom.node
*elements
*set
ring.geom.elem.platens
*set
ring.geom.elem.cyl
