# $Id: makefile,v 1.1 2002-10-05 00:49:22 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = SV_NeoHookean2D.o \
	SV_NeoHookean3D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = SV_NeoHookean2D.cpp \
	SV_NeoHookean3D.cpp

# headers
HDR = SV_NeoHookean2D.h \
	SV_NeoHookean3D.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories

SUB_DIR = 
subdir_driver: 
#SUB_DIR is empty
#subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
