# $Id: makefile,v 1.6 2001-11-28 22:05:45 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ArgSpecT.o \
	CommandSpecT.o \
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
	ArgSpecT.cpp \
	CommandSpecT.cpp \
	iConsoleBaseT.cpp \
	iConsoleObjectT.cpp \
	iConsoleT.cpp

# headers
HDR = \
	ArgSpecT.h \
	CommandSpecT.h \
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
