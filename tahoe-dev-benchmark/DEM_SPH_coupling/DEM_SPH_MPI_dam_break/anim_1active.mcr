#!MC 900

#This macro will animate zones in a dataset 
#such that one or more zones will remain active throughout
#the animation.  You choose the group of zones you want 
#to remain active by entering the numbers of the
#start and end zones of the group.  If you want only
#one zone to remain active, then the start zone should
#be the same as the end zone.  Your group of zones must 
#include the first zone or the last zone in the zone list.

#Last modified by Chris Idso - Tecplot, Inc 12/30/2003

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

$!PROMPTFORTEXTSTRING |AFNAME| 
    INSTRUCTIONS = "Choose a name for the animation file:"

$!PROMPTFORTEXTSTRING |FIRSTZONE|
   INSTRUCTIONS = "What is the number of the first zone of the group that is always active?"

$!PROMPTFORTEXTSTRING |ENDZONE|
   INSTRUCTIONS = "What is the number of the last zone of the group that is always active?"

   
$!EXPORTSETUP 
   EXPORTFNAME = "|AFNAME|.|ANIMFILE|" 
   EXPORTFORMAT = |EXPFORM|

$!IF |FIRSTZONE| == 1      
  $!VARSET |ZONESTOANIMATE| = (|NUMZONES| - |ENDZONE|)
  $!VARSET |ZONESTOADD| = |ENDZONE|
$!ENDIF   
   
$!IF |ENDZONE| == |NUMZONES| 
   $!VARSET |ZONESTOADD| = 0      
   $!VARSET |ZONESTOANIMATE| = (|FIRSTZONE| - 1)
$!ENDIF
   

$!LOOP |ZONESTOANIMATE|
  $!VARSET |REALZONE| = (|LOOP|+|ZONESTOADD|)
  $!ACTIVEFIELDZONES = [|FIRSTZONE|-|ENDZONE|]
  $!ACTIVEFIELDZONES += [|REALZONE|]    
  $!REDRAW
  $!IF |LOOP| == 1
     $!EXPORTSTART
  $!ENDIF
  $!IF |LOOP| <> 1
    $!EXPORTNEXTFRAME       
  $!ENDIF
  
$!ENDLOOP
$!EXPORTFINISH

$!REMOVEVAR |ZONESTOANIMATE|
$!REMOVEVAR |ZONESTOADD|
$!REMOVEVAR |FIRSTZONE|
$!REMOVEVAR |ENDZONE|
$!REMOVEVAR |REALZONE|
$!REMOVEVAR |ANIMFORMAT|
$!REMOVEVAR |ANIMFILE|
$!REMOVEVAR |EXPFORM|
