# $Id: makefile,v 1.7 2001-04-16 17:30:52 rjones Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MultiplierContact2DT.o \
	MultiplierContact3DT.o \
	ContactElementT.o \
	ContactSearchT.o \
	ContactSurfaceT.o \
	ContactNodeT.o \
	SurfaceT.o \
	FaceT.o \
	LineL2FaceT.o \
	QuadL4FaceT.o \
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
	MultiplierContact2DT.cpp \
	MultiplierContact3DT.cpp \
	ContactElementT.cpp \
	ContactSearchT.cpp \
	ContactSurfaceT.cpp \
	ContactNodeT.cpp \
	SurfaceT.cpp \
	FaceT.cpp \
	LineL2FaceT.cpp \
	QuadL4FaceT.cpp \
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
        QuadL4FaceT.h \
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
