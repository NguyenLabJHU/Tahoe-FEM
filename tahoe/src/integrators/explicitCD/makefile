# $Id: makefile,v 1.6 2003-01-29 09:10:59 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eExplicitCD.o \
	ExplicitCDIntegrator.o \
	nExplicitCD.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eExplicitCD.cpp \
	ExplicitCDIntegrator.cpp \
	nExplicitCD.cpp

# headers
HDR = \
	ExplicitCD.h \
	eExplicitCD.h \
	ExplicitCDIntegrator.h \
	nExplicitCD.h

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
