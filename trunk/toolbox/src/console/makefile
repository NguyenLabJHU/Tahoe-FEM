# $Id: makefile,v 1.5 2001-02-22 23:45:19 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	iConsoleBaseT.o \
	iConsoleObjectT.o \
	iConsoleT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	iConsoleBaseT.cpp \
	iConsoleObjectT.cpp \
	iConsoleT.cpp

# headers
HDR = \
	iConsoleBaseT.h \
	iConsoleObjectT.h \
	iConsoleT.h
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
