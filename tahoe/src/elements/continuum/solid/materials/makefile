# $Id: makefile,v 1.18 2003-01-27 07:00:27 paklein Exp $
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

# materials configuration files
include ../config/MaterialsConfig.make

# subdirectories
SUB_DIR = \
	primitives \
	managers \
	interfaces \
	ABAQUS_UMAT \
	ABAQUS_VUMAT \
	fossum \
	Hookean \
	Ogden \
	plasticity_DP \
	plasticity_J2NlHard \
	rod \
	Simo \
	nlcsolver \
	tevp \
	$(MATERIALS_CONFIG_DIRECTORIES)

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
