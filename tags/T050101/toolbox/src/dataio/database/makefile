# $Id: makefile,v 1.5 2001-02-22 23:45:21 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AbaqusT.o \
	EnSightT.o \
	ExodusT.o \
	ModelFileT.o \
	TecPlotT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	AbaqusT.cpp \
	EnSightT.cpp \
	ExodusT.cpp \
	ModelFileT.cpp \
	TecPlotT.cpp

# headers
HDR = \
	AbaqusT.h \
	EnSightT.h \
	ExodusT.h \
	ModelFileT.h \
	TecPlotT.h
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
