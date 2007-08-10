*version
1.0
*title
simple 2D wave BSM benchmark
*dimensions
96   # number of nodes
2   # number of spatial dimensions
2   # number of element sets
# [ID] [nel] [nen] 1 filled, 2 CB elements
1  25  4		
2  50  4		
4   # number of node sets 1 top, 2 bottom, 3 left, 4 right
1  6
2  6
3  16
4  16
# [ID] [nnd]
0   # number of side sets
# end dimensions
*nodesets
*set	# top layer of nodes
./FE_geom/topnodes.dat
*set	# bottom layer of nodes
./FE_geom/bottomnodes.dat
*set	# left layer of nodes
./FE_geom/left.dat
*set	# right layer of nodes
./FE_geom/right.dat
*sidesets
*elements
*set
./FE_geom/ntfmesh.geom.es0
*set
./FE_geom/cbmesh.geom.es0
*nodes
./FE_geom/femesh.geom.nd
