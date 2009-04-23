# $Id: makefile,v 1.1 2009-04-23 03:03:52 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = 	SSOptimize_MatT.o \
	SSOptimize_MatList3DT.o \
	SSKStV_Optimize.o \
	SSVisco_Optimize.o \
	SSVE_Opto_test.o \
	FSFiberOptimize_MatT.o \
	FSFiberOptimize_MatListT.o \
	AnisoCorneaVisco_Opt.o 


OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = 	SSOptimize_MatT.cpp \
	SSOptimize_MatList3DT.cpp \
	SSKStV_Optimize.cpp \
	SSVisco_Optimize.cpp \
	SSVE_Opto_test.cpp \
	FSFiberOptimize_MatT.cpp \
	FSFiberOptimize_MatListT.cpp \
	AnisoCorneaVisco_Opt.cpp  

# headers
HDR = 	SSOptimize_MatT.h \
	SSOptimize_MatList3DT.h \
	SSKStV_Optimize.h \
	SSVisco_Optimize.h \
	SSVE_Opto_test.h \
	FSFiberOptimize_MatT.h \
	FSFiberOptimize_MatListT.h \
	AnisoCorneaVisco_Opt.h  

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
#SUB_DIR = 
#subdir_driver: subdir_loop

# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
