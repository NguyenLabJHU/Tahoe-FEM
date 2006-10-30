# $Id: makefile,v 1.1 2006-10-30 23:32:05 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
#OBJ = 
#OBJ_LINK = $(OBJ)
#DEP = $(OBJ:.o=.d)
OBJ = dummy
OBJ_LINK =
DEP = $(MACRO_DIR)/empty

# sources
#SRC = 

# headers
#HDR = 
#HDR_LINK = $(HDR:.h=.h_link)

# element configuration file
include $(SRC_DIR)/config/DevelopmentMaterialsConfig.make

# subdirectories
SUB_DIR = \
	Ogden \
	Reese-Govindjee \
	Simo-Holzapfel \
	SMP

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
