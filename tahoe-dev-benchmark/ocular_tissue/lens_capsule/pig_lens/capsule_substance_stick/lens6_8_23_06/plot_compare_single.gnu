#set term post color
#set term post eps color 22
#set output "plotpt1mod.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
plot 'datapt1.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 3, 'Rpt1.txt'  t "finite element model" with lines 1
pause -1 "Click on this window and hit return when you are done"
#plot 'datapt3.txt' using 1:2 t "experimental data for 0.3 mm/s" with lines 3, 'Rpt3.txt'  t "finite element model" with lines 1
#pause -1 "Click on this window and hit return when you are done"
#plot 'datapt5.txt' using 1:2 t "experimental data for 0.5 mm/s" with lines 1, 'Rpt5.txt'  t "finite element model" with lines 3
#pause -1 "Click on this window and hit return when you are done"
#
