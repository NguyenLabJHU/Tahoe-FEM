# $Id: makefile,v 1.4 2001-02-22 23:46:18 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ContinuumElementT.o \
	DiffusionT.o \
	ElasticT.o \
	LocalizerT.o \
	TotLag_FDElasticT.o \
	UpLag_FDElasticT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ContinuumElementT.cpp \
	DiffusionT.cpp \
	ElasticT.cpp \
	LocalizerT.cpp \
	TotLag_FDElasticT.cpp \
	UpLag_FDElasticT.cpp

# headers
HDR = \
	ContinuumElementT.h \
	DiffusionT.h \
	ElasticT.h \
	LocalizerT.h \
	TotLag_FDElasticT.h \
	UpLag_FDElasticT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	meshfree
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
