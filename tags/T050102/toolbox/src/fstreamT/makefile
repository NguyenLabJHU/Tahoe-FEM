# $Id: makefile,v 1.6 2002-01-05 06:55:01 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ifstreamT.o \
	ofstreamT.o \
	fstreamT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ifstreamT.cpp \
	ofstreamT.cpp \
	fstreamT.cpp

# headers
HDR = \
	fstreamT.h \
	ifstreamT.h \
	ofstreamT.h
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
