#!/bin/csh

# $argv[1] is .in from Dakota
# $argv[2] is .out returned to Dakota

# echo $PATH

#---------------
# configuration 
#---------------
#set TAHOE         = mpirun -np 2 tahoe
#set TAHOE          = ${HOME}/Tahoe/tahoe/tahoe
set TAHOE          = /home/soils/cven/faculty/regueiro/Code/Tahoe/serial/tahoe-install/tahoe/tahoe_s
#set DAKOTA_UTIL    = ${HOME}/Tahoe/contrib/Dakota
set DAKOTA_UTIL    = /home/soils/cven/faculty/regueiro/Code/Tahoe/serial/tahoe-install/contrib/Dakota
set EXTRACT        = extract_1D
set TAHOE_TEMPLATE = lens_compress.xml.tmpl
set RESULT_EXT     = io4.exo
set ERROR_FUNC     = least_square
set REF_RESULTS    = 2_10_03_Eye02.txt

#--------------
# create tahoe input file
#--------------

echo "* creating input file"
perl ${DAKOTA_UTIL}/dakota2tahoe.pl $argv[1] ${TAHOE_TEMPLATE} $argv[1].xml

#-----
# run 
#-----

echo "* running analysis"
${TAHOE} -f $argv[1].xml > $argv[1].console

#---------------
# extract data
#---------------

echo "* extracting data"
${DAKOTA_UTIL}/${EXTRACT} $argv[1].${RESULT_EXT} $argv[1].${RESULT_EXT} $argv[1].fvst_Y

#---------------
# compute error function
#---------------

echo "* computing objective function"
${DAKOTA_UTIL}/${ERROR_FUNC} ${REF_RESULTS} $argv[1].fvst_Y $argv[2]

#---------------------------------------------------------------
# clean up
#---------------------------------------------------------------
# comment out these lines if running dakota in parallel

#echo "* cleaning up"
#rm -f "$argv[1]".io*
#rm -f "$argv[1]".out "$argv[1]".console
#rm -f "$argv[1]".echo.xml "$argv[1]".valid.xml
#rm -f $argv[1].*xml 
#rm -f $argv[1].fvst_Y
