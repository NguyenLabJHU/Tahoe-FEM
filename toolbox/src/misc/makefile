# $Id: makefile,v 1.7 2002-01-09 18:25:25 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FileCrawlerT.o \
	GroupAverageT.o \
	LocalParabolaT.o \
	MonitorT.o \
	Rotate2DT.o \
	SecantMethodT.o \
	SpectralDecompT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FileCrawlerT.cpp \
	GroupAverageT.cpp \
	LocalParabolaT.cpp \
	MonitorT.cpp \
	Rotate2DT.cpp \
	SecantMethodT.cpp \
	SpectralDecompT.cpp

# headers
HDR = \
	FileCrawlerT.h \
	GroupAverageT.h \
	LocalParabolaT.h \
	MonitorT.h \
	Rotate2DT.h \
	SecantMethodT.h \
	SpectralDecompT.h \
	MapSetT.h

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
