# $Id: makefile,v 1.6 2002-04-02 23:19:24 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eStaticController.o \
	nStaticController.o \
	nLinearStaticController.o \
	LinearStaticController.o \
	StaticController.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eStaticController.cpp \
	nStaticController.cpp \
	nLinearStaticController.cpp \
	LinearStaticController.cpp \
	StaticController.cpp

# headers
HDR = \
	StaticT.h \
	eStaticController.h \
	nStaticController.h \
	nLinearStaticController.h \
	LinearStaticController.h \
	StaticController.h
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
