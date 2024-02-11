#set term post color
#set term post eps color 22
#set output "stress.eps"
#set term post eps
#set title "STRESS STRAIN curves"
#set output "react_uni_bw.eps"
#set xrange [0:.3]
#set yrange [0:1.4e5]
#set key .15,1.5e8
#set xlabel "gamma"
#set ylabel "stress"
#plot 'out.dat' using 3:2  t "ODE solution" with lines 1, 'eldata_mesh1_1.txt' using 2:5  t "1 element FE solution" with lines 2, 'eldata_mesh10_1.txt' using 2:5  t "100 element FE solution" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
set key 700,1.5e8
set xlabel "time"
set ylabel "stress"
plot 'out.dat' using 1:2  t "ODE solution" with lines 1, 'eldata_mesh1_1.txt' using 1:5  t "1 element FE solution" with lines 2, 'eldata_mesh10_1.txt' using 1:5  t "100 element FE solution" with lines 3, 's_alpha.txt' using 1:2  t "100 vector element FE solution" with lines 4
pause -1 "Click on this window and hit return when you are done"
#
#set xlabel "STRAIN"
#set ylabel "GAMMA"
#plot 'out.dat' using 3:4  t "APS model" with lines 1
#pause -1 "Click on this window and hit return when you are done"
#
#set yrange [0:9e4]
#set key .1,8.4e7
#set xlabel "gamma"
#set ylabel "KAPPA (Pa)"
#set ylabel "kappa"
#plot 'out.dat' using 3:5  t "ODE solution" with lines 1, 'eldata_mesh1_1.txt' using 2:8  t "1 element FE solution" with lines 2, 'eldata_mesh10_1.txt' using 2:8  t "100 element FE solution" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
#set xrange [0:1000]
#set key 500,1.5e5
#set xlabel "displ"
#set ylabel "reaction"
#plot 'R_mesh1.txt' using 1:2  t "1 element regular mesh" with lines 1, 'R_mesh4irreg.txt' using 1:2  t "4 element irregular mesh" with lines 4, 'R_mesh_irreg.txt' using 1:2  t "36 element irregular mesh" with lines 3, 'R_mesh10.txt' using 1:2  t "100 element regular mesh" with lines 2
#pause -1 "Click on this window and hit return when you are done"
#
