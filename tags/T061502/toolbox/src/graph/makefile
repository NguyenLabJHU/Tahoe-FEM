# $Id: makefile,v 1.5 2001-02-22 23:45:24 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CMReLabellerT.o \
	GraphBaseT.o \
	GraphT.o \
	MSRBuilderT.o \
	PartitionT.o \
	PriorityQueueT.o \
	ReLabellerT.o \
	RootedLevelT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CMReLabellerT.cpp \
	GraphBaseT.cpp \
	GraphT.cpp \
	MSRBuilderT.cpp \
	PartitionT.cpp \
	PriorityQueueT.cpp \
	ReLabellerT.cpp \
	RootedLevelT.cpp

# headers
HDR = \
	CMReLabellerT.h \
	GraphBaseT.h \
	GraphT.h \
	MSRBuilderT.h \
	PartitionT.h \
	PriorityQueueT.h \
	ReLabellerT.h \
	RootedLevelT.h
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
