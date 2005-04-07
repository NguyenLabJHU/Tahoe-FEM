#!/bin/csh

# $argv[1] is .in from Dakota
# $argv[2] is .out returned to Dakota

# echo $PATH

#---------------
# configuration 
#---------------
set TAHOE          = tahoe

set DAKOTA_UTIL    = ../
set EXTRACT        = extract_1D
set TAHOE_TEMPLATE = cornea.xml.tmpl
set RESULT_EXT     = io1.run
set ERROR_FUNC     = least_square
set REF_RESULTS    = experimental_data

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
${DAKOTA_UTIL}/${EXTRACT} $argv[1].${RESULT_EXT} $argv[1].${RESULT_EXT} $argv[1].fvst_X

#---------------
# compute error function
#---------------

echo "* computing objective function"
${DAKOTA_UTIL}/${ERROR_FUNC} ${REF_RESULTS} $argv[1].fvst_X $argv[2]

#---------------------------------------------------------------
# clean up
#---------------------------------------------------------------

#echo "* cleaning up"
#rm -f $argv[1].*xml $argv[1].io* $argv[1].out $argv[1].console $argv[1].fvst_X
