# $Id: makefile,v 1.16 2003-03-02 18:58:31 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ACME_Contact3DT.o \
	AugLagContact2DT.o \
	AugLagContact3DT.o \
	ContactT.o \
	Contact2DT.o \
	Contact3DT.o \
	Contact3DT_K.o \
	PenaltyContact2DT.o \
	PenaltyContact3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ACME_Contact3DT.cpp \
	AugLagContact2DT.cpp \
	AugLagContact3DT.cpp \
	ContactT.cpp \
	Contact2DT.cpp \
	Contact3DT.cpp \
	Contact3DT_K.cpp \
	PenaltyContact2DT.cpp \
	PenaltyContact3DT.cpp

# headers
HDR = \
	ACME_Contact3DT.h \
	AugLagContact2DT.h \
	AugLagContact3DT.h \
	ContactT.h \
	Contact2DT.h \
	Contact3DT.h \
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
