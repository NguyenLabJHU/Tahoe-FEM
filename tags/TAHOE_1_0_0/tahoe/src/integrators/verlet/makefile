# $Id: makefile,v 1.2 2003-01-29 09:11:03 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eVerlet.o \
	VerletIntegrator.o \
	nVerlet.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eVerlet.cpp \
	VerletIntegrator.cpp \
	nVerlet.cpp

# headers
HDR = \
	Verlet.h \
	eVerlet.h \
	VerletIntegrator.h \
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
