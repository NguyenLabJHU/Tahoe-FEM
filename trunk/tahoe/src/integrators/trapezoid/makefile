# $Id: makefile,v 1.5 2002-04-02 23:19:25 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eTrapezoid.o \
	nTrapezoid.o \
	TrapezoidController.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eTrapezoid.cpp \
	nTrapezoid.cpp \
	TrapezoidController.cpp

# headers
HDR = \
	eTrapezoid.h \
	nTrapezoid.h \
	Trapezoid.h \
	TrapezoidController.h

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
