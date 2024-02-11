# plotting instructions for gnuplot
set title 'check of F = 1 + 0.1 t'
set yrange[0.009:0.01]
set xlabel 'time'
set ylabel 'l(1,1)'
plot 'velocity_gradient.dat' using 1:2 with points title 'tahoe', \
'velocity_gradient.dat' using 1:(0.01)/(1+0.01*$1) with lines title 'exact'
