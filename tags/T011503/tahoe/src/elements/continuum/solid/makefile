# $Id: makefile,v 1.12 2002-11-22 23:34:07 creigh Exp $
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
	SimoFiniteStrainT.o \
	SimoQ1P0.o \
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
	SimoFiniteStrainT.cpp \
	SimoQ1P0.cpp \
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
	SimoFiniteStrainT.h \
	SimoQ1P0.h \
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
