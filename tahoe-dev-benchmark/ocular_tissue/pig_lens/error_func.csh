#!/bin/csh 

# $argv[1] is .in from Dakota
# $argv[2] is .out returned to Dakota

# run two experimental comparisons concurrently
set COMP1     = ../pig_lens2a
set COMP2     = ../pig_lens3

echo "* COMPARISON $COMP1 "
cp -f  $argv[1]  $COMP1
(cd $COMP1 ; ./error_func.csh $argv[1] $argv[2])


echo "* COMPARISON $COMP2 "
cp -f  $argv[1]  $COMP2
(cd $COMP2 ; ./error_func.csh $argv[1] $argv[2])

echo "* SUMMING ERROR "
if (-e ./error_add) then
./error_add $COMP1/$argv[2] $COMP2/$argv[2] 
else
  echo "perl script error_add needs to exist and be executable\n";
endif

mv error_sum $argv[2]
