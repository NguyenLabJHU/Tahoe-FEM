# $Id: makefile,v 1.1 2003-01-29 09:10:46 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ContinuumElementT.o \
	ContinuumMaterialT.o \
	MaterialListT.o \
	MaterialSupportT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ContinuumElementT.cpp \
	ContinuumMaterialT.cpp \
	MaterialListT.cpp \
	MaterialSupportT.cpp

# headers
HDR = \
	ContinuumElementT.h \
	ContinuumMaterialT.h \
	MaterialListT.h \
	MaterialSupportT.h

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
