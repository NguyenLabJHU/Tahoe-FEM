# $Id: makefile,v 1.5 2003-01-30 00:43:38 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	EAM.o \
	EAMFCC2D.o \
	EAMFCC3D.o \
	EAMFCC3DMatT.o \
	EAMFCC3DSym.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	EAM.cpp \
	EAMFCC2D.cpp \
	EAMFCC3D.cpp \
	EAMFCC3DMatT.cpp \
	EAMFCC3DSym.cpp

# headers
HDR = \
	EAMFCC2D.h \
	EAMFCC3D.h \
	EAMFCC3DMatT.h \
	EAMFCC3DSym.h \
	EAM.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	EAM_glue
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
