# $Id: makefile,v 1.4 2001-02-22 23:46:35 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IsoVIB2D.o \
	J2IsoVIB2DLinHardT.o \
	OgdenIsoVIB2D.o \
	VIB2D.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	IsoVIB2D.cpp \
	J2IsoVIB2DLinHardT.cpp \
	OgdenIsoVIB2D.cpp \
	VIB2D.cpp

# headers
HDR = \
	IsoVIB2D.h \
	J2IsoVIB2DLinHardT.h \
	OgdenIsoVIB2D.h \
	VIB2D.h
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
