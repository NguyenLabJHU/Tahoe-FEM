# $Id: makefile,v 1.4 2005-08-12 05:27:11 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GPart.basics.o \
	GPart.DDsepInfo.o \
	GPart.DDviaFishnet.o \
	GPart.DDviaProjection.o \
	GPart.domSegMap.o \
	GPart.identifyWideSep.o \
	GPart.init.o \
	GPart.makeYCmap.o \
	GPart.RBviaDDsep.o \
	GPart.smoothBisector.o \
	GPart.smoothBy2layers.o \
	GPart.smoothYSep.o \
	GPart.split.o \
	GPart.TwoSetViaBKL.o \
	GPart.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	GPart.basics.c \
	GPart.DDsepInfo.c \
	GPart.DDviaFishnet.c \
	GPart.DDviaProjection.c \
	GPart.domSegMap.c \
	GPart.identifyWideSep.c \
	GPart.init.c \
	GPart.makeYCmap.c \
	GPart.RBviaDDsep.c \
	GPart.smoothBisector.c \
	GPart.smoothBy2layers.c \
	GPart.smoothYSep.c \
	GPart.split.c \
	GPart.TwoSetViaBKL.c \
	GPart.util.c

# headers
HDR = \
	DDsepInfo.h \
	GPart.h
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
