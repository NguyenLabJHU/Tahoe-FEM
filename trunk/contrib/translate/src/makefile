# $Id: makefile,v 1.4 2001-09-21 15:49:25 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	TranslateIOManager.o \
	ExtractIOManager.o \
	PointPlots.o \
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
	PointPlots.cpp \
	main.cpp

# headers
HDR = \
	TranslateIOManager.h \
	ExtractIOManager.h \
	PointPlots.h \
	DataManagerT.h

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
