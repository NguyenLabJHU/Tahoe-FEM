# $Id: makefile,v 1.6 2001-07-19 13:56:17 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AbaqusT.o \
	EnSightT.o \
	ExodusT.o \
	ModelFileT.o \
	TecPlotT.o \
	PatranT.o
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
	TecPlotT.cpp \
	PatranT.cpp

# headers
HDR = \
	AbaqusT.h \
	EnSightT.h \
	ExodusT.h \
	ModelFileT.h \
	TecPlotT.h \
	PatranT.h
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
