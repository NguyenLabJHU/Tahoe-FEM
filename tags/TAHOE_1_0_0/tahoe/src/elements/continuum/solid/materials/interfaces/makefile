# $Id: makefile,v 1.3 2003-01-29 09:10:51 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SolidMatSupportT.o \
	FSMatSupportT.o \
	SSMatSupportT.o \
	D2FSMatSupportT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	SolidMatSupportT.cpp \
	FSMatSupportT.cpp \
	SSMatSupportT.cpp \
	D2FSMatSupportT.cpp

# headers
HDR = \
	SolidMatSupportT.h \
	FSMatSupportT.h \
	SSMatSupportT.h \
	D2FSMatSupportT.h

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
