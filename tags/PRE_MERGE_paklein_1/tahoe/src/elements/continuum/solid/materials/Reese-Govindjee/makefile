# $Id: makefile,v 1.1 2002-10-04 23:59:00 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = RGBaseT.o \
	RG_VDSplit2D.o \
	RG_VDSplit3D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = RGBaseT.cpp \
	RG_VDSplit2D.cpp \
	RG_VDSplit3D.cpp

# headers
HDR = RGBaseT.h \
	RG_VDSplit2D.h \
	RG_VDSplit3D.h


HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = VIB \
	Mooney-Rivlin
subdir_driver: subdir_loop

# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
