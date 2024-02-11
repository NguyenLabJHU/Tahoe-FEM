#set term post color
set term post eps color 22
set output "plotpt1.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
#
#plot 'datapt1_w_precond.txt' using 2:1 t "data for 0.1 mm/s, showing preconditioning" with lines 1
plot 'datapt1_w_precond.txt' using 2:1 not with lines 3
pause -1 "Click on this window and hit return when you are done"
#
#plot 'datapt1.txt' using 1:2 t "experimental data for 0.1 mm/s without preconditioning in plot" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
#plot 'datapt3_w_precond.txt' using 2:1 t "data for 0.3 mm/s, showing preconditioning" with lines 1
#plot 'datapt3_w_precond.txt' using 2:1 not with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
#plot 'datapt3.txt' using 1:2 t "experimental data for 0.3 mm/s without preconditioning in plot" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
#plot 'datapt5_w_precond.txt' using 2:1 t "data for 0.5 mm/s, showing preconditioning" with lines 1
#plot 'datapt5_w_precond.txt' using 2:1 not with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
#plot 'datapt5.txt' using 1:2 t "experimental data for 0.5 mm/s without preconditioning in plot" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
