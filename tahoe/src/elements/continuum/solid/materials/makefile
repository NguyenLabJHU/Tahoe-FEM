# $Id: makefile,v 1.21 2003-01-30 00:43:37 paklein Exp $
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
include $(SRC_DIR)/config/SolidMaterialsConfig.make

# subdirectories
SUB_DIR = \
	primitives \
	managers \
	interfaces \
	Hookean \
	Ogden \
	Simo \
	nlcsolver \
	$(SOLID_MATERIALS_CONFIG_DIRECTORIES)

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
