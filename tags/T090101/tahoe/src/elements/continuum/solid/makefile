# $Id: makefile,v 1.8 2001-08-20 15:43:13 paklein Exp $
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
	NonsingularContinuumT.o \
	SimoFiniteStrainT.o \
	SmallStrainT.o \
	TotalLagrangianT.o \
	UpdatedLagrangianT.o \
	UpLagr_ExternalFieldT.o

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
	NonsingularContinuumT.cpp \
	SimoFiniteStrainT.cpp \
	SmallStrainT.cpp \
	TotalLagrangianT.cpp \
	UpdatedLagrangianT.cpp \
	UpLagr_ExternalFieldT.cpp

# headers
HDR = \
	ContinuumElementT.h \
	DiffusionT.h \
	ElasticT.h \
	FiniteStrainT.h \
	LocalizerT.h \
	NonsingularContinuumT.h \
	SimoFiniteStrainT.h \
	SmallStrainT.h \
	TotalLagrangianT.h \
	UpdatedLagrangianT.h \
	UpLagr_ExternalFieldT.h

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
