# $Id: makefile,v 1.2.2.1 2002-12-19 03:09:14 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CommunicatorT.o \
	MessageT.o \
	AllGatherT.o \
	PointToPointT.o	

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CommunicatorT.cpp \
	MessageT.cpp \
	AllGatherT.cpp \
	PointToPointT.cpp

# headers
HDR = \
	CommunicatorT.h \
	MessageT.h \
	AllGatherT.h \
	PointToPointT.h

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
