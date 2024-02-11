#set term post color
#set term post eps color 22
#set title "STRESS STRAIN curves"
#set output "sigeps.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "STRAIN"
set ylabel "STRESS"
plot 'eldataloc_1_2pt5.txt' using ($9):($3) t "2.5/sec" with lines 1, 'eldataloc_1_pt25.txt' using 9:3  t "0.25/sec" with lines 2, 'eldataloc_1_pt025.txt' using 9:3  t "0.025/sec" with lines 3
pause -1 "Click on this window and hit return when you are done"
#
