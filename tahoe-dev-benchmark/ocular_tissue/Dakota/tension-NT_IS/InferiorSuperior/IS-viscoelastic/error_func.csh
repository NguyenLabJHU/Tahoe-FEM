#!/bin/csh

# $argv[1] is .in from Dakota
# $argv[2] is .out returned to Dakota

# echo $PATH

#---------------
# configuration 
#---------------
set TAHOE          = tahoe

set DAKOTA_UTIL    = ${TAHOE_HOME}/contrib/Dakota
set EXTRACT        = extract_1D
set TAHOE_TEMPLATE = tension-strip.xml.tmpl
set RESULT_EXT     = io1.run
set ERROR_FUNC     = least_square
set REF_RESULTS    = tension-reaction.dat

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

echo "* cleaning up"
rm -f $argv[1].console 
#rm -f $argv[1].out 
rm -f $argv[1].io* 
#rm -f $argv[1].*xml
#rm -f $argv[1].fvst_X 
cp -f $argv[1].fvst_X current.fvst_X
cp -f $argv[1].xml current.xml
cp -f $argv[1] current.in
#cp $argv[1] $argv[1].parameters
#cp $argv[2] $argv[2].norm
