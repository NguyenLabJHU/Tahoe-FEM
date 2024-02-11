#set term post color
#set term post eps color 22
#set term post eps
#set output "react_uni_bw.eps"
#set xrange [0:.3]
#set yrange [0:1.4e5]
#set key .15,1.5e8
#
set xlabel "D_Y"
set ylabel "T_n"
plot 'node_rigid_5.txt' using 3:7  t "T_n vs D_Y" with lines 1
pause -1 "Click on this window and hit return when you are done"
#
set xlabel "d_n"
plot 'node_rigid_5.txt' using 5:7  t "T_n vs d_n" with lines 1
pause -1 "Click on this window and hit return when you are done"
#
set xlabel "up_n"
plot 'node_rigid_5.txt' using 9:7  t "T_n vs up_n" with lines 1
pause -1 "Click on this window and hit return when you are done"
#
