*version
1.0
*title
MD only BSM 2D LJ hex benchmark nearest neighbor
*dimensions
218   # number of nodes 
2   # number of spatial dimensions
1   # number of element sets
# [ID] [nel] [nen]
1	218	1
7   # number of node sets 
# [ID] [nnd]
1	14	# bottom
2	14	# top
3	8	# left
4	8	# right
5	15	# ghost bottom
6	15	# ghost top
7	188	# real
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
*set	# "ghost atoms"
./MD_geom/ghost_bottom.dat
*set	# "ghost atoms"
./MD_geom/ghost_top.dat
*set	# real
./MD_geom/real.dat
*sidesets
*elements
*set
./MD_geom/allatoms.es0
*nodes	# all atom (including ghost) initial coordinates given
./MD_geom/allatoms.nd