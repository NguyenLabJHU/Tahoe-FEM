# $Id: makefile,v 1.16 2003-03-26 20:04:24 cjkimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InelasticDuctile2DT.o \
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
	MR2DT.cpp \
	MR_RP2DT.cpp

# headers
HDR = \
	InelasticDuctile2DT.h \
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
