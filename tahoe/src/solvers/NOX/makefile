# $Id: makefile,v 1.2 2002-07-03 23:11:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	NOX_Tahoe_Vector.o \
	NOX_Tahoe_Group.o \
	NOXSolverT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	NOXSolverT.cpp \
	NOX_Tahoe_Group.cpp \
	NOX_Tahoe_Vector.cpp

# headers
HDR = \
	NOXSolverT.h \
	NOX_Tahoe_Group.h \
	NOX_Tahoe_Vector.h

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
