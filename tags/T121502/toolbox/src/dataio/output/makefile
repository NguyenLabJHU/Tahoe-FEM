# $Id: makefile,v 1.8 2002-08-01 16:33:34 saubry Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AVSOutputT.o \
	AbaqusOutputT.o \
	EnSightOutputT.o \
	ExodusOutputT.o \
	FE_ASCIIT.o \
	OutputBaseT.o \
	OutputSetT.o \
	TecPlotOutputT.o \
	PatranOutputT.o \
	ParaDynOutputT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	AVSOutputT.cpp \
	AbaqusOutputT.cpp \
	EnSightOutputT.cpp \
	ExodusOutputT.cpp \
	FE_ASCIIT.cpp \
	OutputBaseT.cpp \
	OutputSetT.cpp \
	TecPlotOutputT.cpp \
	PatranOutputT.cpp \
	ParaDynOutputT.cpp
# headers
HDR = \
	AVSOutputT.h \
	AbaqusOutputT.h \
	EnSightOutputT.h \
	ExodusOutputT.h \
	FE_ASCIIT.h \
	OutputBaseT.h \
	OutputSetT.h \
	TecPlotOutputT.h \
	PatranOutputT.h \
	ParaDynOutputT.h
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
