#set term post color
#set term post eps color 22
#set term post eps
#set title "STRESS STRAIN curves"
#set output "react_uni_bw.eps"
#set xrange [0:.3]
#set yrange [0:1.4e5]
#set key .15,1.5e8
set xlabel "gamma"
set ylabel "stress"
plot 'node_3.txt' using 3:2  t "solution" with lines 1
#plot 'node_3.txt' using 3:2  t "solution" with lines 1, 'node_3_change.txt' using 3:2  t "solution" with lines 1
pause -1 "Click on this window and hit return when you are done"
#
