# $Id: makefile,v 1.1 2002-10-04 23:59:02 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = RGVIB2D.o \
	ViscVIB.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = RGVIB2D.cpp \
	ViscVIB.cpp

# headers
HDR = RGVIB2D.h \
	ViscVIB.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = ViscosityFunctions
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
