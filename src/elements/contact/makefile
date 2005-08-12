# $Id: makefile,v 1.19 2005-08-12 05:37:58 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ContactElementT.o \
	ContactFaceT.o \
	ContactLineL2FaceT.o \
	ContactLineQ3FaceT.o \
	ContactNodeT.o \
	ContactQuadL4FaceT.o \
	ContactSearchT.o \
	ContactSurfaceT.o \
	FaceT.o \
	LineL2FaceT.o \
	LineQ3FaceT.o \
	MultiplierContactElement3DT.o \
	MultiplierContactElement2DT.o \
	FrictionalContactElement2DT.o \
	PenaltyContactElement2DT.o \
	PenaltyContactElement3DT.o \
	QuadL4FaceT.o \
	SurfaceT.o \
	TriaL3FaceT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	ContactElementT.cpp \
	ContactFaceT.cpp \
	ContactLineL2FaceT.cpp \
	ContactLineQ3FaceT.cpp \
	ContactNodeT.cpp \
	ContactQuadL4FaceT.cpp \
	ContactSearchT.cpp \
	ContactSurfaceT.cpp \
	FaceT.cpp \
	LineL2FaceT.cpp \
	LineQ3FaceT.cpp \
	MultiplierContactElement3DT.cpp \
	MultiplierContactElement2DT.cpp \
	FrictionalContactElement2DT.cpp \
	PenaltyContactElement2DT.cpp \
	PenaltyContactElement3DT.cpp \
	QuadL4FaceT.cpp \
	SurfaceT.cpp \
	TriaL3FaceT.cpp

# headers
HDR = \
	ContactElementT.h \
	ContactFaceT.h \
	ContactLineL2FaceT.h \
	ContactLineQ3FaceT.h \
	ContactNodeT.h \
	ContactQuadL4FaceT.h \
	ContactSearchT.h \
	ContactSurfaceT.h \
	FaceT.h \
	LineL2FaceT.h \
	LineQ3FaceT.h \
	MultiplierContactElement3DT.h \
	MultiplierContactElement2DT.h \
	FrictionalContactElement2DT.h \
	PenaltyContactElement2DT.h \
	PenaltyContactElement3DT.h \
	QuadL4FaceT.h \
	SurfaceT.h \
	TriaL3FaceT.h \
	vector2D.h \
	vector3D.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
