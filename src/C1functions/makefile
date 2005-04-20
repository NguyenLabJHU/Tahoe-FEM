# $Id: makefile,v 1.27 2005-04-20 23:47:27 thao Exp $
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
#DEP = /dev/null

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
DEPEND = /dev/null
include $(DEPEND)
