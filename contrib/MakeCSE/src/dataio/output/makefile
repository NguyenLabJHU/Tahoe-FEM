# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ExodusOutputT.o \
	FE_ASCIIT.o \
	OutputBaseT.o \
	OutputSetT.o \
	EnSightOutputT.o \
	AbaqusOutputT.o \
	TecPlotOutputT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ExodusOutputT.cpp \
	FE_ASCIIT.cpp \
	OutputBaseT.cpp \
	OutputSetT.cpp \
	EnSightOutputT.cpp \
	AbaqusOutputT.cpp \
	TecPlotOutputT.cpp

# headers
HDR = \
	ExodusOutputT.h \
	FE_ASCIIT.h \
	OutputBaseT.h \
	OutputSetT.h \
	EnSightOutputT.h \
	AbaqusOutputT.h \
	TecPlotOutputT.h

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
