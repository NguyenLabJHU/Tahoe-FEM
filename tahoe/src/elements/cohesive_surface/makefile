# $Id: makefile,v 1.4 2001-02-22 23:46:16 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CSEAnisoT.o \
	CSEBaseT.o \
	CSEIsoT.o \
	MeshFreeCSEAnisoT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CSEAnisoT.cpp \
	CSEBaseT.cpp \
	CSEIsoT.cpp \
	MeshFreeCSEAnisoT.cpp

# headers
HDR = \
	CSEAnisoT.h \
	CSEBaseT.h \
	CSEIsoT.h \
	MeshFreeCSEAnisoT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	surface_potentials
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
