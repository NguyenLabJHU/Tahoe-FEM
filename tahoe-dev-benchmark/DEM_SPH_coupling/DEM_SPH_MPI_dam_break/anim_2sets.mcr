#!MC 900
#Last modified by Chris Idso, Tecplot, Inc  12/30/2003
#This macro is intended to address the case where
#there are two sets of zones, each with n timesteps,
#and they are both loaded to the same frame in Tecplot.
#This macro can animate this data such that the first image
#of the animation contains the first zone of set 1 and the 
#first zone of set 2, the second image contains the second
#zone of set 1and the second zone of set 2, etc.

#For example:  Two sets, each with 10 timesteps and nothing
#else loaded.
#The first animation image would contain zones 1 and 11
#the first zone of each set, the second image would contain
#zones 2 and 12, the second zone of each set, then
#zones 3 and 13 and so on, until
#zones 10 and 20, the last image.

#This macro will prompt for the number of the first zone
#of both sets,, and the number of timesteps.

$!PROMPTFORTEXTSTRING |SETONEZONE|
    INSTRUCTIONS = "Enter the number of the first zone of set 1:"

$!PROMPTFORTEXTSTRING |SETTWOZONE|
    INSTRUCTIONS = "Enter the number of the first zone of set 2:"

$!PROMPTFORTEXTSTRING |NUMINSET|
    INSTRUCTIONS = "Enter the number of zones in a set:"

$!PROMPTFORTEXTSTRING |AFNAME|
    INSTRUCTIONS = "Choose a name for the animation file:"

$!VARSET |ANIMFILE| = "swf"
$!VARSET |EXPFORM| = "FLASH"

$!PROMPTFORTEXTSTRING |ANIMFORMAT|
   INSTRUCTIONS = "Enter animation format (1 - FLASH, 2 - AVI, 3 - RM):"

$!IF |ANIMFORMAT| == 2
  $!VARSET |ANIMFILE| = "avi"
  $!VARSET |EXPFORM| = "AVI"
$!ENDIF 

$!IF |ANIMFORMAT| == 3
  $!VARSET |ANIMFILE| = "rm"
  $!VARSET |EXPFORM| = "RASTERMETAFILE"
$!ENDIF 

$!EXPORTSETUP 
  EXPORTFNAME = "|AFNAME|.|ANIMFILE|" 
  EXPORTFORMAT = |EXPFORM|

$!LOOP |NUMINSET|
  $!ACTIVEFIELDZONES = [|SETONEZONE|,|SETTWOZONE|]
  $!REDRAW
  $!IF |LOOP| == 1
    $!EXPORTSTART      
  $!ENDIF
  $!IF |LOOP| <> 1
    $!EXPORTNEXTFRAME 
  $!ENDIF
  $!VARSET |SETONEZONE| += 1
  $!VARSET |SETTWOZONE| += 1
$!ENDLOOP

$!EXPORTFINISH

$!REMOVEVAR |SETONEZONE|
$!REMOVEVAR |SETTWOZONE|
$!REMOVEVAR |NUMINSET|
$!REMOVEVAR |EXPFORM|
$!REMOVEVAR |ANIMFORMAT|
$!REMOVEVAR |ANIMFILE|
