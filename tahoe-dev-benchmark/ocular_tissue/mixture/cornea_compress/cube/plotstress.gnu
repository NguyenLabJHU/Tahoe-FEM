#set term post color
#set term post eps color 22
#set title "STRESS STRAIN curves"
#set output "sigeps.eps"
#
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
set xlabel "TIME (sec)"
set ylabel "STRESS (Pa)"
plot 'ndata_1.txt' using 1:abs(3)  t "D = 3.06e-12 sec" with lines 1, 'ndata_1_x10e5.txt' using 1:3  t "D = 3.06e-5 sec" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
set xlabel "TIME (sec)"
set ylabel "C_F_ref (kg/m^3)"
plot 'ndata_1.txt' using 1:4  t "D = 3.06e-12 sec" with lines 1, 'ndata_1_x10e5.txt' using 1:4  t "D = 3.06e-5 sec" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
