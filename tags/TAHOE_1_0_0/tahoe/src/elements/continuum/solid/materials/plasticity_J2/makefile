# $Id: makefile,v 1.8 2002-04-17 23:58:18 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	J2PrimitiveT.o \
	J2QLLinHardT.o \
	J2QL2DLinHardT.o \
	J2Simo2D.o \
	J2Simo3D.o \
	J2SSKStV.o \
	J2SSKStV2D.o \
	J2SSC0HardeningT.o \
	J2SimoC0HardeningT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	J2PrimitiveT.cpp \
	J2QLLinHardT.cpp \
	J2QL2DLinHardT.cpp \
	J2Simo2D.cpp \
	J2Simo3D.cpp \
	J2SSKStV.cpp \
	J2SSKStV2D.cpp \
	J2SSC0HardeningT.cpp \
	J2SimoC0HardeningT.cpp

# headers
HDR = \
	J2PrimitiveT.h \
	J2QLLinHardT.h \
	J2QL2DLinHardT.h \
	J2Simo2D.h \
	J2Simo3D.h \
	J2SSKStV.h \
	J2SSKStV2D.h \
	J2SSC0HardeningT.h \
	J2SimoC0HardeningT.h

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
