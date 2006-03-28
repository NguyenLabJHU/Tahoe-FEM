# $Id: makefile,v 1.19 2006-03-28 16:50:49 raregue Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InelasticDuctile2DT.o \
	InelasticDuctile_RP2DT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	InelasticDuctile2DT.cpp \
	InelasticDuctile_RP2DT.cpp

# headers
HDR = \
	InelasticDuctile2DT.h \
	InelasticDuctile_RP2DT.h

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
