#set term post color
set term post eps color 22
#set term post eps
set output "R.eps"
#set xrange [0:.01]
set yrange [0:2.25]
#set key .15,1.5e8
#set xlabel "up_n"
#set ylabel "T_n"
set xlabel "displacement"
set ylabel "force"
#plot 'node_5.txt' using 3:5  t "solution" with lines 1
plot 'R.txt' using 1:2 not with lines 1
pause -1 "Click on this window and hit return when you are done"
#
