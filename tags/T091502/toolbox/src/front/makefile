# $Id: makefile,v 1.5 2001-02-22 23:45:23 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FrontNodeT.o \
	FrontSegmentT.o \
	FrontT.o \
	Front2DT.o \
	Front3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FrontNodeT.cpp \
	FrontSegmentT.cpp \
	FrontT.cpp \
	Front2DT.cpp \
	Front3DT.cpp

# headers
HDR = \
	FrontNodeT.h \
	FrontSegmentT.h \
	FrontT.h \
	Front2DT.h \
	Front3DT.h
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
