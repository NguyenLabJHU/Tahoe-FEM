# $Id: makefile,v 1.6 2001-05-17 18:05:37 ebmarin Exp $
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
	CB \
	continuum \
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
	plasticity_macro \
	nlcsolver \
	tevp

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
