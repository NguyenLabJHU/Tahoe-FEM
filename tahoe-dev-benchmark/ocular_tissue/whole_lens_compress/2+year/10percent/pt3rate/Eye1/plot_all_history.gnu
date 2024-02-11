#set term post color
#set term post eps color 22
#set output "plot.eps"
#set xrange [0:.1]
#set yrange [0:1.4e5]
#set key .4,1e5
#
set xlabel "TIME (sec)"
set ylabel "FORCE (N)"
#plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'R.txt'  t "model" with lines 2
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.1.fvst_Y'  t "model1" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.2.fvst_Y'  t "model2" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.3.fvst_Y'  t "model3" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.4.fvst_Y'  t "model4" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.5.fvst_Y'  t "model5" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.6.fvst_Y'  t "model6" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.7.fvst_Y'  t "model7" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.8.fvst_Y'  t "model8" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.9.fvst_Y'  t "model9" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.10.fvst_Y'  t "model10" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.11.fvst_Y'  t "model11" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.12.fvst_Y'  t "model12" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.13.fvst_Y'  t "model13" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.14.fvst_Y'  t "model14" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.15.fvst_Y'  t "model15" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.16.fvst_Y'  t "model16" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.17.fvst_Y'  t "model17" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.18.fvst_Y'  t "model18" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.19.fvst_Y'  t "model19" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.20.fvst_Y'  t "model20" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.21.fvst_Y'  t "model21" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.22.fvst_Y'  t "model22" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.23.fvst_Y'  t "model23" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.24.fvst_Y'  t "model24" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.25.fvst_Y'  t "model25" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.26.fvst_Y'  t "model26" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
plot 'Eye1Low.txt' using 1:2 t "experimental data for 0.1 mm/s" with lines 1, 'error_func.in.27.fvst_Y'  t "model27" with lines 2
pause -1 "Click on this window and hit return when you are done"
#
