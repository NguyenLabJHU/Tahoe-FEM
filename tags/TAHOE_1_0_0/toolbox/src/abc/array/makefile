# $Id: makefile,v 1.6 2002-10-08 21:52:38 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ArraySettings.o \
	dArrayT.o \
	dArray2DT.o \
	iArrayT.o \
	iArray2DT.o \
	iAutoArrayT.o \
	sArrayT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ArraySettings.cpp \
	dArrayT.cpp \
	dArray2DT.cpp \
	iArrayT.cpp \
	iArray2DT.cpp \
	iAutoArrayT.cpp \
	sArrayT.cpp
# headers
HDR = \
	ArrayT.h \
	Array2DT.h \
	AutoArrayT.h \
	dArrayT.h \
	dArray2DT.h \
	iArrayT.h \
	iArray2DT.h \
	iAutoArrayT.h \
	nArrayT.h \
	nArray2DT.h \
	sArrayT.h \
	RaggedArray2DT.h \
	VariRaggedArray2DT.h
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
