# $Id: makefile,v 1.4 2001-02-22 23:46:09 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	eLinearHHTalpha.o \
	eNLHHTalpha.o \
	HHTalpha.o \
	LinearHHTalpha.o \
	NLHHTalpha.o \
	nLinearHHTalpha.o \
	nNLHHTalpha.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	eLinearHHTalpha.cpp \
	eNLHHTalpha.cpp \
	HHTalpha.cpp \
	LinearHHTalpha.cpp \
	NLHHTalpha.cpp \
	nLinearHHTalpha.cpp \
	nNLHHTalpha.cpp

# headers
HDR = \
	eLinearHHTalpha.h \
	eNLHHTalpha.h \
	HHTalpha.h \
	LinearHHTalpha.h \
	NLHHTalpha.h \
	nLinearHHTalpha.h \
	nNLHHTalpha.h
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
