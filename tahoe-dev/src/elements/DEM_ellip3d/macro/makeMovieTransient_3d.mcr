#!MC 1410
$!VarSet |MFBD| = ''
$!Varset |NAME| = "3d"
$!OPENLAYOUT  "|NAME|.lay"

$!ExportSetup ExportFormat = PNG
$!ExportSetup ImageWidth = 2000
$!ExportSetup ExportFName = '|NAME|_ini.png'
$!Export 
  ExportRegion = CurrentFrame
  
$!ExportSetup ExportFormat = MPEG4
$!ExportSetup ExportFName = '|NAME|.mp4'
$!AnimateTime 
  StartTime = 1
  EndTime = 100
  Skip = 1
  CreateMovieFile = Yes

$!ExportSetup ExportFormat = PNG
$!ExportSetup ImageWidth = 2000
$!ExportSetup ExportFName = '|NAME|_end.png'
$!Export 
  ExportRegion = CurrentFrame
