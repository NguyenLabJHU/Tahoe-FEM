#!MC 750

##
## This macro will copy a zone a specified
## number of times, then rotate each copy around a point.
## Great for making propellers, turbines, flowers, 
## pinwheels etc.   
## This macro leaves the plot in 2d mode.  
## 10/21/99 Dave Smith
##last edited by Chris Idso 8/26/03


$!PromptforTextString |ZONE|
  Instructions = "Which zone do you want to duplicate and rotate"
$!PromptforTextString |Blades|
  Instructions = "How many duplicates do you want made?"
$!Loop |Blades|
  $!DUPLICATEZONE 
    SOURCEZONE = |ZONE|
  $!Varset |I| += 1
  $!RenameDataSetZone
    Zone = |NumZones|
    Name = "Blade #|Loop|"
$!EndLoop
$!ActiveFieldZones -= [|ZONE|]
$!Varset |FirstZone| = (|NumZones| - |Blades| + 1)
$!Varset |I| = |FirstZone| 
$!Varset |RotationAngle| = (360/|Blades|)
$!FRAMEMODE = TWOD
$!While |I| < |NumZones|
  $!Rotate2d
    Angle = |RotationAngle|
    ZoneList = [|FirstZone|-|I|]
$!Varset |I| += 1
$!EndWhile
$!VIEW DATAFIT
