# $Id: makefile,v 1.7 2003-01-29 09:11:01 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eStaticIntegrator.o \
	nStaticIntegrator.o \
	nLinearStaticIntegrator.o \
	LinearStaticIntegrator.o \
	StaticIntegrator.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eStaticIntegrator.cpp \
	nStaticIntegrator.cpp \
	nLinearStaticIntegrator.cpp \
	LinearStaticIntegrator.cpp \
	StaticIntegrator.cpp

# headers
HDR = \
	StaticT.h \
	eStaticIntegrator.h \
	nStaticIntegrator.h \
	nLinearStaticIntegrator.h \
	LinearStaticIntegrator.h \
	StaticIntegrator.h

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
