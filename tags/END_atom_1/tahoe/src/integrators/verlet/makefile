# $Id: makefile,v 1.1 2002-06-27 23:55:57 hspark Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eVerlet.o \
	VerletController.o \
	nVerlet.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eVerlet.cpp \
	VerletController.cpp \
	nVerlet.cpp

# headers
HDR = \
	Verlet.h \
	eVerlet.h \
	VerletController.h \
	nVerlet.h
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
