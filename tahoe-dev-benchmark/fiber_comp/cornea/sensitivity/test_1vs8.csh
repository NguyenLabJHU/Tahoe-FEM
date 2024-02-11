tahoe -f rate350_1.xml
tahoe -f rate350_8.xml
rm -f 1.uvst_Y 8.uvst_Y
extract_1D rate350_1.io1.exo rate350_1.io1.exo 1.uvst_Y
extract_1D rate350_8.io1.exo rate350_8.io1.exo 8.uvst_Y

echo 'plot "1.uvst_Y" w lp, "8.uvst_Y" w l; pause 4' | gnuplot -
