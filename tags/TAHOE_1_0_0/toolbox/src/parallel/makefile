# $Id: makefile,v 1.3 2003-01-27 06:42:48 paklein Exp $
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
