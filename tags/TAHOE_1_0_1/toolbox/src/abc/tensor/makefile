# $Id: makefile,v 1.2 2001-10-15 17:20:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	dTensor3DT.o \
	dTensor4DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	dTensor3DT.cpp \
	dTensor4DT.cpp

# headers
HDR = \
	TensorT.h \
	Tensor3DT.h \
	Tensor4DT.h \
	dTensor3DT.h \
	dTensor4DT.h

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
