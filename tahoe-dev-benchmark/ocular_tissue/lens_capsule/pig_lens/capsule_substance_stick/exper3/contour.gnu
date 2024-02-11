#
#set term post
#set output "contour.ps"
set term post eps 16
#
# Set misc parameters
set cntrparam cubicspline
#set cntrparam bspline
set cntrparam points 5
#set cntrparam order 3
#set cntrparam levels discrete 1,1.5,2,2.5
#set cntrparam levels auto 5
#set cntrparam levels incremental 0,.05,.4
# delete next two lines to see surface, else is at base
#set contour base
#set nosurface
set noclip
set parametric
set data style lines
set view 45, 150
#
# Set title and labels
#set title "excavation"
set noclabel
set hidden3d
#set xlabel "X1, m"
#set ylabel "X2, m" -3,0
#set tics out
set ticslevel 0
set noxtics
set noytics
#set xtics 0,.025,.15
#set ytics 0,.025,.05
#set xrange[-30:30]
#set yrange[25:35]
#
# Plot the data file
#set output "contour.eps"
splot 'params1.out' not
pause -1 "Click on this window and hit return when you are done"
#
#set zrange[-50:150]
#set zrange[-150:200]
#set zrange[0:50]
#set zrange[0:0.5]
#set zlabel "p" ,-1
#set zlabel "tauoct" ,-1
#set zlabel "gammoct" ,-1
#set zlabel "effpl" ,-1
#set output "excc9crsemidp.eps"
#splot 'params1.out' not
#pause -1 "Click on this window and hit return when you are done"
#
#set zrange[0:1.5]
#set zrange[0:.15]
#set zlabel "GAMMAOCT" ,-1
#set output "ccrseenhgam.eps"
#splot 'ccrseenhgam.out' not
#pause -1 "Click on this window and hit return when you are done"
#set output "cfineenhgam.eps"
#splot 'cfineenhgam.out' not
#pause -1 "Click on this window and hit return when you are done"
#set output "cfinerenhgam.eps"
#splot 'cfinerenhgam.out' not
#pause -1 "Click on this window and hit return when you are done"
#set zrange[0:.4]
#set zrange[0:2.5]
#set zlabel "EFFPLSTRN" ,-1
#set output "ccrseenheff.eps"
#splot 'ccrseenheff.out' not
#pause -1 "Click on this window and hit return when you are done"
#set output "cfineenheff.eps"
#splot 'cfineenheff.out' not
#pause -1 "Click on this window and hit return when you are done"
#set output "cfinerenheff.eps"
#splot 'cfinerenheff.out' not
#pause -1 "Click on this window and hit return when you are done"
#
# To try to mix two contour plots like for the excavation problem
#set term table
#set nosurface
#set output 'contour.dat'
#set term x11
#plot "contour.dat"
#
