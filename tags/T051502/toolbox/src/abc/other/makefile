# $Id: makefile,v 1.9 2002-03-19 18:59:20 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	dRangeArrayT.o \
	LocalArrayT.o \
	StringT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	dRangeArrayT.cpp \
	LocalArrayT.cpp \
	StringT.cpp

# headers
HDR = \
	AutoFill2DT.h \
	dRangeArrayT.h \
	LocalArrayT.h \
	pArrayT.h \
	pMatrixT.h \
	RowAutoFill2DT.h \
	StringT.h \
	Vector3T.h \
	Vector2T.h
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
