# $Id: makefile,v 1.6 2003-01-29 09:11:02 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eTrapezoid.o \
	nTrapezoid.o \
	TrapezoidIntegrator.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eTrapezoid.cpp \
	nTrapezoid.cpp \
	TrapezoidIntegrator.cpp

# headers
HDR = \
	eTrapezoid.h \
	nTrapezoid.h \
	Trapezoid.h \
	TrapezoidIntegrator.h

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
