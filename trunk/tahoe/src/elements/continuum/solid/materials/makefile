# $Id: makefile,v 1.10 2001-07-18 21:29:43 paklein Exp $
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
	ABAQUS_UMAT \
	ABAQUS_VUMAT \
	CB \
	EAM \
	Hookean \
	managers \
	modCBSW \
	Ogden \
	plasticity_DP \
	plasticity_J2 \
	primitives \
	rod \
	Simo \
	VIB \
	plasticity_crystl \
	plasticity_macro \
	nlcsolver \
	tevp \
	SimoVisco \
	VISCVIB

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
