# $Id: makefile,v 1.6 2003-01-29 09:10:58 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IntegratorT.o \
	eIntegratorT.o \
	nIntegratorT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	IntegratorT.cpp \
	eIntegratorT.cpp \
	nIntegratorT.cpp

# headers
HDR = \
	IntegratorT.h \
	eIntegratorT.h \
	nIntegratorT.h

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
