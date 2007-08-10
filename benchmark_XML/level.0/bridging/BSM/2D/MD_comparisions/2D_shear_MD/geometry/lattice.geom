*version
1.0
*title
MD 2D LJ hex benchmark nearest neighbor
*dimensions
536   # number of nodes 
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1	536	1
5   # number of node sets 
# [ID] [nnd]
1	14	# bottom
2	14	# top
3	18	# left
4	18	# right
5	536	# real
0   # number of side sets
# end dimensions
*nodesets
*set	# bottom atoms
./MD_geom/bottom.dat
*set	# top atoms
./MD_geom/top.dat
*set	#left
./MD_geom/left.dat
*set	#right
./MD_geom/right.dat
*set	# real
./MD_geom/real.dat
*sidesets
*elements
*set
./MD_geom/allatoms.es0
*nodes	# all atom (including ghost) initial coordinates given
./MD_geom/allatoms.nd
