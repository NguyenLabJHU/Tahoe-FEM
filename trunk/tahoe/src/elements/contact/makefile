# $Id: makefile,v 1.11 2002-01-23 00:01:53 rjones Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	PenaltyContactElement2DT.o \
	AdhesionContact2DT.o \
	MultiplierContact2DT.o \
	MultiplierContact3DT.o \
	ContactElementT.o \
	ContactSearchT.o \
	ContactSurfaceT.o \
	ContactNodeT.o \
	SurfaceT.o \
	FaceT.o \
	LineL2FaceT.o \
	LineQ3FaceT.o \
	QuadL4FaceT.o \
	ContactFaceT.o \
	ContactLineL2FaceT.o \
	ContactLineQ3FaceT.o \
	ContactQuadL4FaceT.o \
	ACME_Contact3DT.o \
	AugLagContact2DT.o \
	ContactT.o \
	Contact2DT.o \
	Contact3DT.o \
	DOFElementT.o \
	PenaltyContact2DT.o \
	PenaltyContact3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	PenaltyContactElement2DT.cpp \
	AdhesionContact2DT.cpp \
	MultiplierContact2DT.cpp \
	MultiplierContact3DT.cpp \
	ContactElementT.cpp \
	ContactSearchT.cpp \
	ContactSurfaceT.cpp \
	ContactNodeT.cpp \
	SurfaceT.cpp \
	FaceT.cpp \
	LineL2FaceT.cpp \
	LineQ3FaceT.cpp \
	QuadL4FaceT.cpp \
	ContactFaceT.cpp \
	ContactLineL2FaceT.cpp \
	ContactLineQ3FaceT.cpp \
	ContactQuadL4FaceT.cpp \
	ACME_Contact3DT.cpp \
	AugLagContact2DT.cpp \
	ContactT.cpp \
	Contact2DT.cpp \
	Contact3DT.cpp \
	DOFElementT.cpp \
	PenaltyContact2DT.cpp \
	PenaltyContact3DT.cpp

# headers
HDR = \
	PenaltyContactElement2DT.h \
	AdhesionContact2DT.h \
	MultiplierContact2DT.h \
	MultiplierContact3DT.h \
	ContactSurfaceT.h \
	ContactElementT.h \
	ContactSearchT.h \
	ContactSurfaceT.h \
	ContactNodeT.h \
	SurfaceT.h \
	FaceT.h \
	LineL2FaceT.h \
	LineQ3FaceT.h \
	QuadL4FaceT.h \
	ContactFaceT.cpp \
	ContactLineL2FaceT.cpp \
	ContactLineQ3FaceT.cpp \
	ContactQuadL4FaceT.cpp \
	vector2D.h \
	vector3D.h \
	ACME_Contact3DT.h \
	AugLagContact2DT.h \
	ContactT.h \
	Contact2DT.h \
	Contact3DT.h \
	DOFElementT.h \
	PenaltyContact2DT.h \
	PenaltyContact3DT.h
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
