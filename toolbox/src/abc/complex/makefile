# $Id: makefile,v 1.6 2001-12-19 01:41:45 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ComplexT.o \
	zArrayT.o \
	zMatrixT.o \
	zTensor3DT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ComplexT.cpp \
	zArrayT.cpp \
	zMatrixT.cpp \
	zTensor3DT.cpp

# headers
HDR = \
	ComplexT.h \
	zArrayT.h \
	zMatrixT.h \
	zTensor3DT.h

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
