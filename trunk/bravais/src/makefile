# $Id: makefile,v 1.7 2003-06-06 23:16:03 saubry Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	PeriodicTableT.o \
	PerTabEntryT.o \
	VolumeT.o \
	CrystalLatticeT.o \
	BoxT.o \
	AsperityT.o \
	CUBT.o \
	FCCT.o \
	BCCT.o \
	DIAT.o \
	HEXT.o \
	OutPutLatticeT.o \
	MeshAtom.o \
	MakeCrystalT.o

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
	VolumeT.cpp \
        CrystalLatticeT.cpp \
        BoxT.cpp \
	AsperityT.cpp \
	CUBT.cpp \
	FCCT.cpp \
	BCCT.cpp \
	DIAT.cpp \
	HEXT.cpp \
	OutPutLatticeT.cpp \
	MeshAtom.cpp \
	MakeCrystalT.cpp
 
# headers
HDR = \
	PeriodicTableT.h \
	PerTabEntryT.h \
	VolumeT.h \
        CrystalLatticeT.h \
        BoxT.h \
	AsperityT.h \
	CUBT.h \
	FCCT.h \
	BCCT.h \
	DIAT.h \
	HEXT.h \
	OutPutLatticeT.h \
	MeshAtom.h \
	MakeCrystalT.h

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
