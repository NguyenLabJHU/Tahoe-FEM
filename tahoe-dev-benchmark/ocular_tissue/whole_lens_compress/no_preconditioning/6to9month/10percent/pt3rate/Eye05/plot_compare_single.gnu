#set term post eps color 22
#set output "plotpt3modpt1.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
plot '6-9_10_03_Eye05.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 3, 'R.txt'  t "finite element model" with lines 1
pause -1 "Click on this window and hit return when you are done"
