# $Id: makefile,v 1.8 2002-06-08 20:20:53 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ElementCardT.o \
	ElementMatrixT.o \
	FBC_CardT.o \
	IC_CardT.o \
	KBC_CardT.o \
	ScheduleT.o \
	SWDataT.o \
	TimeSequence.o \
	Traction_CardT.o \
	ElementBlockDataT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ElementCardT.cpp \
	ElementMatrixT.cpp \
	FBC_CardT.cpp \
	IC_CardT.cpp \
	KBC_CardT.cpp \
	ScheduleT.cpp \
	SWDataT.cpp \
	TimeSequence.cpp \
	Traction_CardT.cpp \
	ElementBlockDataT.cpp

# headers
HDR = \
	ElementCardT.h \
	ElementMatrixT.h \
	FBC_CardT.h \
	IC_CardT.h \
	KBC_CardT.h \
	ScheduleT.h \
	SWDataT.h \
	TimeSequence.h \
	Traction_CardT.h \
	ElementBlockDataT.h 

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	globalmatrix
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
