# $Id: makefile,v 1.3 2002-10-05 03:44:33 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = SSSV_KStV2D.o \
	SSSV_KStV3D.o \
	FDSV_KStV2D.o \
	FDSV_KStV3D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = SSSV_KStV2D.cpp \
	SSSV_KStV3D.cpp \
	FDSV_KStV2D.cpp \
	FDSV_KStV3D.cpp

# headers
HDR = SSSV_KStV2D.h \
	SSSV_KStV3D.h \
	FDSV_KStV2D.h \
	FDSV_KStV3D.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR =
# subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
