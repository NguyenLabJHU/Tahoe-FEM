# $Id: makefile,v 1.13 2003-01-29 09:10:49 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SolidElementT.o \
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
	SolidElementT.cpp \
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
	SolidElementT.h \
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
	materials \
	meshfree 

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
