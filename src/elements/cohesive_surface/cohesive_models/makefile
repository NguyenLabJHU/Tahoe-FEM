# $Id: makefile,v 1.17 2003-08-08 16:03:18 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InelasticDuctile2DT.o \
	InelasticDuctile_RP2DT.o \
	MR2DT.o \
	MR_RP2DT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	InelasticDuctile2DT.cpp \
	InelasticDuctile_RP2DT.cpp \
	MR2DT.cpp \
	MR_RP2DT.cpp

# headers
HDR = \
	InelasticDuctile2DT.h \
	InelasticDuctile_RP2DT.h \
	MR2DT.h \
	MR_RP2DT.h

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
