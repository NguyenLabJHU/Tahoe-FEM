# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AbaqusInputT.o \
	EnSightInputT.o \
	ExodusInputT.o \
	InputBaseT.o \
	TextInput.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	AbaqusInputT.cpp \
	EnSightInputT.cpp \
	ExodusInputT.cpp \
	InputBaseT.cpp \
	TextInput.cpp

# headers
HDR = \
	AbaqusInputT.h \
	EnSightInputT.h \
	ExodusInputT.h \
	InputBaseT.h \
	TextInput.h

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
