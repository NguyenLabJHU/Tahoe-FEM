# $Id: makefile,v 1.9 2002-08-01 16:35:05 saubry Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AVST.o \
	AbaqusResultsT.o \
	EnSightT.o \
	ExodusT.o \
	ModelFileT.o \
	TecPlotT.o \
	PatranT.o \
	ParaDynT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	AVST.cpp \
	AbaqusResultsT.cpp \
	EnSightT.cpp \
	ExodusT.cpp \
	ModelFileT.cpp \
	TecPlotT.cpp \
	PatranT.cpp \
	ParaDynT.cpp

# headers
HDR = \
	AVST.h \
	AbaqusResultsT.h \
	AbaqusVariablesT.h \
	EnSightT.h \
	ExodusT.h \
	ModelFileT.h \
	TecPlotT.h \
	PatranT.h \
	ParaDynT.h
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
