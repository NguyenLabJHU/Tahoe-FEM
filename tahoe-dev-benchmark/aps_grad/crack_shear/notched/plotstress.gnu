#set term post color
#set term post eps color 22
#set title "STRESS STRAIN curves"
#set output "sigeps.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "TIME (sec)"
set ylabel "REACTION (N)"
#plot 'R.txt' using 1:2  t "no curl" with lines 1, 'R_curl.txt' using 1:2  t "with curl" with lines 2, 'R_curl_zero.txt' using 1:2  t "with curl zero" with lines 3, 'R_zero.txt' using 1:2  t "with curl zero" with lines 4
plot 'R_nocurl1.txt' using 1:2  t "no curl - mesh1" with lines 1, 'R_curl1.txt' using 1:2  t "with curl - mesh1" with lines 2, 'R_nocurl2.txt' using 1:2  t "no curl - mesh2" with lines 3, 'R_curl2.txt' using 1:2  t "with curl - mesh2" with lines 4
pause -1 "Click on this window and hit return when you are done"
plot 'R_nocurl1.txt' using 1:2  t "no curl - mesh1" with lines 1, 'R_curl1.txt' using 1:2  t "with curl - mesh1" with lines 2, 'R_nocurl2.txt' using 1:2  t "no curl - mesh2" with lines 3, 'R_curl2.txt' using 1:2  t "with curl - mesh2" with lines 4, 'R_nocurl3.txt' using 1:2  t "no curl - mesh3" with lines 5, 'R_curl3.txt' using 1:2  t "with curl - mesh3" with lines 6
pause -1 "Click on this window and hit return when you are done"
#
