# $Id: makefile,v 1.1 2004-03-20 02:48:12 raregue Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MRPrimitiveT.o \
	MRSSKStV.o \
	MRSSKStV2D.o \
	MRSSNLHardT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MRPrimitiveT.cpp \
	MRSSKStV.cpp \
	MRSSKStV2D.cpp \
	MRSSNLHardT.cpp

# headers
HDR = \
	MRPrimitiveT.h \
	MRSSKStV.h \
	MRSSKStV2D.h \
	MRSSNLHardT.h

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
