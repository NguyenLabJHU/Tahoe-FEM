# $Id: makefile,v 1.5 2001-07-03 01:35:21 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IsoVIB3D.o \
	J2IsoVIB3DLinHardT.o \
	OgdenIsoVIB3D.o \
	VIB3D.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	IsoVIB3D.cpp \
	J2IsoVIB3DLinHardT.cpp \
	OgdenIsoVIB3D.cpp \
	VIB3D.cpp

# headers
HDR = \
	IsoVIB3D.h \
	J2IsoVIB3DLinHardT.h \
	OgdenIsoVIB3D.h \
	VIB3D.h
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
