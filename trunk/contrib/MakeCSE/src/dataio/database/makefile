# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ExodusT.o \
	ModelFileT.o \
	EnSightT.o \
	AbaqusT.o \
	TecPlotT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ExodusT.cpp \
	ModelFileT.cpp \
	EnSightT.cpp \
	AbaqusT.cpp \
	TecPlotT.cpp

# headers
HDR = \
	ExodusT.h \
	ModelFileT.h \
	EnSightT.h \
	AbaqusT.h \
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
