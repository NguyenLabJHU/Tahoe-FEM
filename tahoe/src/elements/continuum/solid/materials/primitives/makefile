# $Id: makefile,v 1.5 2001-07-03 01:35:43 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Anisotropic2DT.o \
	ContinuumMaterialT.o \
	DetCheckT.o \
	DiffusionMaterialT.o \
	FDStructMatT.o \
	IsotropicT.o \
	Material2DT.o \
	NL_E_MatT.o \
	NL_E_Mat2DT.o \
	NL_E_RotMat2DT.o \
	SSStructMatT.o \
	StructuralMaterialT.o \
	ThermalDilatationT.o \
	TensorTransformT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	Anisotropic2DT.cpp \
	ContinuumMaterialT.cpp \
	DetCheckT.cpp \
	DiffusionMaterialT.cpp \
	FDStructMatT.cpp \
	IsotropicT.cpp \
	Material2DT.cpp \
	NL_E_MatT.cpp \
	NL_E_Mat2DT.cpp \
	NL_E_RotMat2DT.cpp \
	SSStructMatT.cpp \
	StructuralMaterialT.cpp \
	ThermalDilatationT.cpp \
	TensorTransformT.cpp

# headers
HDR = \
	Anisotropic2DT.h \
	ContinuumMaterialT.h \
	DetCheckT.h \
	DiffusionMaterialT.h \
	FDStructMatT.h \
	IsotropicT.h \
	Material2DT.h \
	NL_E_MatT.h \
	NL_E_Mat2DT.h \
	NL_E_RotMat2DT.h \
	SSStructMatT.h \
	StructuralMaterialT.h \
	ThermalDilatationT.h \
	TensorTransformT.h

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
