# $Id: makefile,v 1.2 2005-08-12 05:40:38 paklein Exp $
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
#DEP = $(MACRO_DIR)/empty

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
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
