# $Id: makefile,v 1.6 2001-08-03 19:16:43 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InputBaseT.o \
	TahoeInputT.o \
	PatranInputT.o \
	ExodusInputT.o \
	EnSightInputT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	InputBaseT.cpp \
	TahoeInputT.cpp \
	PatranInputT.cpp \
	ExodusInputT.cpp \
	EnSightInputT.cpp
# headers
HDR = \
	InputBaseT.h \
	TahoeInputT.h \
	PatranInputT.h \
	ExodusInputT.h \
	EnSightInputT.h

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
