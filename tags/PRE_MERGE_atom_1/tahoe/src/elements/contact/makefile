# $Id: makefile,v 1.13 2002-06-17 17:15:07 rjones Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	PenaltyContactElement2DT.o \
	PenaltyContactElement3DT.o \
	MultiplierContactElement2DT.o \
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
	TriaL3FaceT.o \
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
	PenaltyContactElement3DT.cpp \
	MultiplierContactElement2DT.cpp \
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
	TriaL4FaceT.o \
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
	PenaltyContactElement3DT.h \
	MultiplierContactElement2DT.h \
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
	TriaL3FaceT.h \
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
