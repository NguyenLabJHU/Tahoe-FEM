# $Id: makefile,v 1.17 2002-11-14 17:05:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
#OBJ = 
#OBJ_LINK = $(OBJ)
#DEP = $(OBJ:.o=.d)
OBJ = dummy
OBJ_LINK =
DEP = /dev/null

# sources
#SRC = 

# headers
#HDR = 
#HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	primitives \
	managers \
	interfaces \
	ABAQUS_UMAT \
	ABAQUS_VUMAT \
	CB \
	EAM \
	elasticity_crystl \
	fossum \
	Hookean \
	modCBSW \
	Ogden \
	plasticity_DP \
	plasticity_J2 \
	plasticity_J2NlHard \
	rod \
	Simo \
	VIB \
	plasticity_crystl \
	plasticity_macro \
	nlcsolver \
	tevp \
	Simo-Holzapfel \
	Reese-Govindjee 

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
