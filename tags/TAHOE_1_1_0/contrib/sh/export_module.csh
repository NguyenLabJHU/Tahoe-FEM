#!/bin/csh
# $Id: export_module.csh,v 1.2 2002-04-04 20:11:19 paklein Exp $
# This shell script creates a tar.gz file from a tahoe module
# USAGE EXAMPLE: export_module.csh [path to toolbox] toolbox.Darwin
# will create toolbox.Darwin.tar.gz from the toolbox module. 
#
set sourceDir=$1
set destDir=$2
if (-e $destDir) then
	echo "$destDir already exists"
	exit 1
endif

# check layout of module
if (! -d $sourceDir/inc) then
	echo "Did not find $sourceDir/inc"
	exit 1
endif

if (-e $sourceDir/makefile) then
	set makef=makefile
else if (-e $sourceDir/Makefile) then
	set makef=Makefile
else
	echo "Did not find $sourceDir/makefile or $sourceDir/Makefile"
	exit 1
endif

# get library name from the makefile
set library=`grep -E "TARGET[[:space:]]*=" $sourceDir/$makef | sed 's/TARGET[[:space:]]*=[[:space:]]*\(.*\)/\1/'`

# check for library
if (! -e $sourceDir/lib/lib$library.a) then
	echo "Did not find library file in $sourceDir/lib/lib$library.a"
	exit 1
endif 

echo "Creating $destDir.tar.gz from directory $sourceDir"	

mkdir $destDir

echo "Copying *.h from $sourceDir/inc"
mkdir $destDir/inc
cp $sourceDir/inc/*.h $destDir/inc
chmod a-w $destDir/inc/*.h

echo "Copying lib$library.a from $sourceDir/lib"
mkdir $destDir/lib
cp $sourceDir/lib/lib$library.a $destDir/lib
chmod a-w $destDir/lib/*.a

echo "Creating tar file"
tar cf $destDir.tar $destDir

echo "Compressing tar file"
gzip -v $destDir.tar

echo "Cleaning up"
rm -rf $destDir

echo "Done."
