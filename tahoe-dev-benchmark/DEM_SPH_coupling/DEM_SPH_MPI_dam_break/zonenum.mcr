#!MC 800
##   This macro creates a new variable and then fills that variable with '1' s
##   for the first zone, '2's for the second zone etc.
##   Author:  Chris Idso - Tecplot,Inc.   

$!VarSet |MFBD| = 'D:\Tec8.0-2'


$!VARSET |NEWVARNUM| = (|NUMVARS| + 1)

$!ALTERDATA 
  EQUATION = 'V|NEWVARNUM| = 1' 

$!LOOP |NUMZONES|
$!ALTERDATA [|LOOP|]
   EQUATION = 'V|NEWVARNUM| = |LOOP|'

$!ENDLOOP

$!RemoveVar |MFBD|
