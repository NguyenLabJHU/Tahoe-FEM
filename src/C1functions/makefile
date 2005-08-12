# $Id: makefile,v 1.28 2005-08-12 05:36:32 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MajumdarBhushan.o \
	GreenwoodWilliamson.o \
	GWPlastic.o \
	PolyDistributionT.o \
	ConHypGeom.o \
	WormLikeChain.o \
	WLCwRep.o \
	FungwRep.o \
	FungType.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	MajumdarBhushan.cpp \
	GreenwoodWilliamson.cpp \
	GWPlastic.cpp \
	PolyDistributionT.cpp \
	ConHypGeom.cpp \
	WormLikeChain.cpp \
	WLCwRep.cpp \
	FungwRep.cpp \
	FungType.cpp

# headers
HDR = \
	MajumdarBhushan.h \
	GreenwoodWilliamson.h \
	GWPlastic.h \
	PolyDistributionT.h \
	ConHypGeom.h \
	WormLikeChain.h \
	WLCwRep.h \
	FungwRep.h \
	FungType.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR =  
#subdir_driver: subdir_loop
#SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
