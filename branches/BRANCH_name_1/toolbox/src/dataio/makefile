# $Id: makefile,v 1.7 2002-01-09 12:09:05 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IOBaseT.o \
	JoinOutputT.o \
	ModelManagerT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	IOBaseT.cpp \
	JoinOutputT.cpp \
	ModelManagerT.cpp

# headers
HDR = \
	IOBaseT.h \
	JoinOutputT.h \
	ModelManagerT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	database \
	input \
	output
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
