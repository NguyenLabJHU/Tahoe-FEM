# $Id: makefile,v 1.4 2001-11-20 01:04:04 recampb Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VTKBodyT.o \
	VTKBodyDataT.o \
	VTKConsoleT.o \
	VTKFrameT.o \
	main.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	VTKBodyT.cpp \
	VTKBodyDataT.cpp \
	VTKConsoleT.cpp \
	VTKFrameT.cpp \
	main.cpp

# headers
HDR = \
	VTKBodyT.h \
	VTKBodyDataT.h \
	VTKConsoleT.h \
	VTKFrameT.h

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
