# $Id: makefile,v 1.1 2003-01-29 09:10:48 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DiffusionMatListT.o \
	DiffusionMatSupportT.o \
	DiffusionMaterialT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DiffusionMatListT.cpp \
	DiffusionMatSupportT.cpp \
	DiffusionMaterialT.cpp

# headers
HDR = \
	DiffusionMatListT.h \
	DiffusionMatSupportT.h \
	DiffusionMaterialT.h

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
