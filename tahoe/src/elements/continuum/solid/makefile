# $Id: makefile,v 1.5 2001-07-03 01:34:53 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ContinuumElementT.o \
	DiffusionT.o \
	ElasticT.o \
	FiniteStrainT.o \
	LocalizerT.o \
	SmallStrainT.o \
	TotalLagrangianT.o \
	UpdatedLagrangianT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ContinuumElementT.cpp \
	DiffusionT.cpp \
	ElasticT.cpp \
	FiniteStrainT.cpp \
	LocalizerT.cpp \
	SmallStrainT.cpp \
	TotalLagrangianT.cpp \
	UpdatedLagrangianT.cpp

# headers
HDR = \
	ContinuumElementT.h \
	DiffusionT.h \
	ElasticT.h \
	FiniteStrainT.h \
	LocalizerT.h \
	SmallStrainT.h \
	TotalLagrangianT.h \
	UpdatedLagrangianT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	meshfree

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
