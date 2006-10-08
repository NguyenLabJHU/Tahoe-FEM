# $Id: makefile,v 1.21 2006-10-08 19:22:19 regueiro Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InelasticDuctile2DT.o \
	InelasticDuctile_RP2DT.o \
	MR2DT.o \
        MR_RP2DT.o \
        MR_NodalRP2DT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	InelasticDuctile2DT.cpp \
	InelasticDuctile_RP2DT.cpp \
	MR2DT.cpp \
        MR_RP2DT.cpp \
        MR_NodalRP2DT.cpp


# headers
HDR = \
	InelasticDuctile2DT.h \
	InelasticDuctile_RP2DT.h \
	MR2DT.h \
        MR_RP2DT.h \
        MR_NodalRP2DT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
