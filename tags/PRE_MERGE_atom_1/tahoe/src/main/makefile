# $Id: makefile,v 1.5 2002-12-05 08:31:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ExecutionManagerT.o \
	FEExecutionManagerT.o \
	FEManagerT.o \
	FEManagerT_mpi.o \
	GlobalT.o \
	IOManager.o \
	IOManager_mpi.o \
	main.o \
	TimeManagerT.o \
	CommManagerT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ExecutionManagerT.cpp \
	FEExecutionManagerT.cpp \
	FEManagerT.cpp \
	FEManagerT_mpi.cpp \
	GlobalT.cpp \
	IOManager.cpp \
	IOManager_mpi.cpp \
	main.cpp \
	TimeManagerT.cpp \
	CommManagerT.cpp

# headers
HDR = \
	ExecutionManagerT.h \
	FEExecutionManagerT.h \
	FEManagerT.h \
	FEManagerT_mpi.h \
	GlobalT.h \
	IOManager.h \
	IOManager_mpi.h \
	TimeManagerT.h \
	CommManagerT.h

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
