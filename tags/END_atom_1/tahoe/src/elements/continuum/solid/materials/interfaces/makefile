# $Id: makefile,v 1.2 2002-11-14 17:06:21 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	D2FDMatSupportT.o \
	DiffusionMatSupportT.o \
	FDMatSupportT.o \
	MaterialSupportT.o \
	SSMatSupportT.o \
	StructuralMatSupportT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	D2FDMatSupportT.cpp \
	DiffusionMatSupportT.cpp \
	FDMatSupportT.cpp \
	MaterialSupportT.cpp \
	SSMatSupportT.cpp \
	StructuralMatSupportT.cpp

# headers
HDR = \
	D2FDMatSupportT.h \
	DiffusionMatSupportT.h \
	FDMatSupportT.h \
	MaterialSupportT.h \
	SSMatSupportT.h \
	StructuralMatSupportT.h

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
