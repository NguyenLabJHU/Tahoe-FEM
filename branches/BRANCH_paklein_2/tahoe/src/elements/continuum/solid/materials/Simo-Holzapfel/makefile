# $Id: makefile,v 1.1 2002-10-05 00:49:20 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = SSSimoViscoT.o \
	FDSimoViscoBaseT.o \
	FDSimoVisco2D.o \
	FDSimoVisco3D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = SSSimoViscoT.cpp \
	FDSimoViscoBaseT.cpp \
	FDSimoVisco2D.cpp \
	FDSimoVisco3D.cpp

# headers
HDR = SSSimoViscoT.h \
	FDSimoViscoBaseT.h \
	FDSimoVisco2D.h \
	FDSimoVisco3D.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = KStV \
	Mooney-Rivlin

subdir_driver: subdir_loop

# SUB_DIR is empty
# subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
