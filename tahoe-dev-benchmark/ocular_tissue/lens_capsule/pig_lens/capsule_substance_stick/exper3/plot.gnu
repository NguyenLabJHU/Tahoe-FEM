#set term post color
#set term post eps color 22
#set output "plot.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
#set title "fit to lens3 data"
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
plot 'experimental_data' t "experimental data" with lines 1, 'error_func.in.6.fvst_Y'  t "model" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
