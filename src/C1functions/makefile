# $Id: makefile,v 1.24 2003-04-25 20:01:39 dzeigle Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MajumdarBhushan.o \
	GreenwoodWilliamson.o \
	ConHypGeom.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MajumdarBhushan.cpp \
	GreenwoodWilliamson.cpp \
	ConHypGeom.cpp

# headers
HDR = \
	MajumdarBhushan.h \
	GreenwoodWilliamson.h \
	ConHypGeom.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR =  
#subdir_driver: subdir_loop
#SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
