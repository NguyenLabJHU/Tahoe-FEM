# $Id: makefile,v 1.7 2002-05-29 12:56:01 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	TranslateIOManager.o \
	ExtractIOManager.o \
	ExtractNode.o \
	ExtractQuad.o \
	PointPlots.o \
	MergeResults.o \
	main.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	TranslateIOManager.cpp \
	ExtractIOManager.cpp \
	ExtractNode.cpp \
	ExtractQuad.cpp \
	PointPlots.cpp \
	MergeResults.cpp \
	main.cpp

# headers
HDR = \
	TranslateIOManager.h \
	ExtractIOManager.h \
	ExtractNode.h \
	ExtractQuad.h \
	MergeResults.h \
	PointPlots.h

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
