# $Id: makefile,v 1.8.2.1 2001-11-13 20:56:01 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	TahoeInputT.o \
	PatranInputT.o \
	ExodusInputT.o \
	EnSightInputT.o \
	AbaqusInputT.o \
	InputFEASCIIT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	TahoeInputT.cpp \
	PatranInputT.cpp \
	ExodusInputT.cpp \
	EnSightInputT.cpp \
	AbaqusInputT.cpp \
	InputFEASCIIT.cpp
# headers
HDR = \
	InputBaseT.h \
	TahoeInputT.h \
	PatranInputT.h \
	ExodusInputT.h \
	EnSightInputT.h \
	AbaqusInputT.h \
	InputFEASCIIT.h

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
