# $Id: makefile,v 1.2 2002-03-06 01:55:43 jzimmer Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	VolumeT.o \
	BoxT.o \
	CrystalLatticeT.o \
	FCCT.o \
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
	BoxT.cpp \
	CrystalLatticeT.cpp \
	FCCT.cpp \
	MakeCrystalT.cpp

# headers
HDR = \
	MakeCrystalT.h \
	VolumeT.h \
	BoxT.h \
	CrystalLatticeT.h \
	FCCT.h

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
