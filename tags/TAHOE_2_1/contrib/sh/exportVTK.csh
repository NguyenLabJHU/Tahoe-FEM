#!/bin/csh
# $Id: exportVTK.csh,v 1.1 2002-04-04 05:37:20 paklein Exp $
# This short shell script creates a tar.gz file from a VTK 4.x
# directory. 
# USAGE EXAMPLE: exportVTK.sh VTK vtk_redhat
# will create vtk_redhat.tar.gz from the VTK directory. 
#
set sourceDir=$1
set destDir=$2
if (-e $destDir) then
	echo "$destDir already exists"
	exit 1
endif

echo "Creating $destDir.tar.gz from directory $sourceDir"	

mkdir $destDir
cp $sourceDir/*.h $destDir

mkdir $destDir/bin
echo "Copying *.a from $sourceDir/bin"
cp $sourceDir/bin/*.a $destDir/bin
chmod a-w $destDir/bin/*.a

set subDirs = (Common Filtering Graphics Imaging IO Parallel Rendering)
foreach n ($subDirs)
	mkdir $destDir/$n
	echo "Copying *.h from $sourceDir/$n"
	cp $sourceDir/$n/*.h $destDir/$n
	chmod a-w $destDir/$n/*.h
end

echo "Creating tar file"
tar cf $destDir.tar $destDir

echo "Compressing tar file"
gzip -v $destDir.tar

echo "Cleaning up"
rm -rf $destDir

echo "Done."
