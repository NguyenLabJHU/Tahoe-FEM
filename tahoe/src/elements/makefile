# $Id: makefile,v 1.11.2.1 2002-12-10 17:08:49 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ElementT.o \
	ElementSupportT.o \
	ElementBaseT.o \
	DOFElementT.o \
	ElementListT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ElementT.cpp \
	ElementSupportT.cpp \
	ElementBaseT.cpp \
	DOFElementT.cpp \
	ElementListT.cpp

# headers
HDR = \
	ElementT.h \
	ElementSupportT.h \
	ElementBaseT.h \
	DOFElementT.h \
	ElementListT.h
HDR_LINK = $(HDR:.h=.h_link)

# element configuration file
include ../config/ElementsConfig.make

# subdirectories
SUB_DIR = \
	BEM_element \
	rod \
	adhesion \
	$(ELEMENTS_CONFIG_DIRECTORIES)

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
