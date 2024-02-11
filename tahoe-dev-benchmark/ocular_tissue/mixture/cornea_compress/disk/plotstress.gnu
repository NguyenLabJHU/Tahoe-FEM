#set term post color
set term post eps color 22
#set title "STRESS STRAIN curves"
set output "force.eps"
#
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
plot 'R.txt' using 1:2  t "reaction force (N)" with lines 1
pause -1 "Click on this window and hit return when you are done"
#
