# $Id: makefile,v 1.1 2002-07-05 22:22:06 hspark Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eGear6.o \
	Gear6Controller.o \
	nGear6.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eGear6.cpp \
	Gear6Controller.cpp \
	nGear6.cpp

# headers
HDR = \
	Gear6.h \
	eGear6.h \
	Gear6Controller.h \
	nGear6.h
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
