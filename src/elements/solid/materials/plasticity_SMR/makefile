# $Id: makefile,v 1.1 2006-07-27 16:26:52 regueiro Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SMRPrimitiveT.o \
	SMRSSKStV.o \
	SMRSSKStV2D.o \
	SMRSSNLHardT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	SMRPrimitiveT.cpp \
	SMRSSKStV.cpp \
	SMRSSKStV2D.cpp \
	SMRSSNLHardT.cpp

# headers
HDR = \
	SMRPrimitiveT.h \
	SMRSSKStV.h \
	SMRSSKStV2D.h \
	SMRSSNLHardT.h

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
