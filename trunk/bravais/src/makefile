# $Id: makefile,v 1.4 2002-07-24 01:15:33 saubry Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	PeriodicTableT.o \
	PerTabEntryT.o \
	MakeCrystalT.o \
	VolumeT.o \
	CrystalLatticeT.o \
	BoxT.o \
	FCCT.o \
	OutPutLatticeT.o \
	MeshAtom.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	main.cpp \
	PeriodicTableT.cpp \
	PerTabEntryT.cpp \
	MakeCrystalT.cpp \
	VolumeT.cpp \
        CrystalLatticeT.cpp \
        BoxT.cpp \
	FCCT.cpp \
	OutPutLatticeT.cpp \
	MeshAtom.cpp

# headers
HDR = \
	PeriodicTableT.h \
	PerTabEntryT.h 

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
