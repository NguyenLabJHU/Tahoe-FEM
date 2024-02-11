#!/bin/csh
set TAHOE = tahoe_d
#
# run test
echo "running ${TAHOE}..."
${TAHOE} -f emmi.1.xml > tahoe.console
#
# extract plot data
echo "extracting plot data..."
# num_ip = 1
grep deform\ meas emmi.1.out | awk -f ip1.awk > velocity_gradient.dat
# num_ip = 8
#grep deform\ meas emmi.1.out | awk -f ip8.awk > velocity_gradient.dat
#
# plot l(1,1)
echo "plotting results with gnuplot..."
gnuplot -persist velocity_gradient.plot
