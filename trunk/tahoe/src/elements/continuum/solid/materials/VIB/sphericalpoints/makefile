# $Id: makefile,v 1.4 2001-02-22 23:46:41 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FCCPtsT.o \
	IcosahedralPtsT.o \
	LatLongPtsT.o \
	SpherePointsT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FCCPtsT.cpp \
	IcosahedralPtsT.cpp \
	LatLongPtsT.cpp \
	SpherePointsT.cpp

# headers
HDR = \
	FCCPtsT.h \
	IcosahedralPtsT.h \
	LatLongPtsT.h \
	SpherePointsT.h
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
