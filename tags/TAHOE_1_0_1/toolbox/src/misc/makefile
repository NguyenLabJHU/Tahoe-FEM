# $Id: makefile,v 1.9 2002-11-25 07:08:30 paklein Exp $
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
	Rotate3DT.o \
	SecantMethodT.o \
	SpectralDecompT.o \
	InverseMapT.o

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
	Rotate3DT.cpp \
	SecantMethodT.cpp \
	SpectralDecompT.cpp \
	InverseMapT.cpp

# headers
HDR = \
	FileCrawlerT.h \
	GroupAverageT.h \
	LocalParabolaT.h \
	MonitorT.h \
	Rotate2DT.h \
	Rotate3DT.h \
	SecantMethodT.h \
	SpectralDecompT.h \
	MapSetT.h \
	InverseMapT.h

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
