# $Id: makefile,v 1.1 2002-10-04 23:59:01 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = RG_NeoHookean2D.o \
	RG_NeoHookean3D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = RG_NeoHookean3D.cpp \
	RG_NeoHookean3D.cpp
# headers
HDR = RG_NeoHookean2D.h \
	RG_NeoHookean3D.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
#SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
