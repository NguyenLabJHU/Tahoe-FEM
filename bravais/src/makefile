# $Id: makefile,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	VolumeT.o \
	LatticeT.o \
	MakeCrystalT.o 
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	main.cpp \
	VolumeT.cpp \
	LatticeT.cpp \
	MakeCrystalT.cpp

# headers
HDR = \
	MakeCrystalT.h \
	VolumeT.h \
	LatticeT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
