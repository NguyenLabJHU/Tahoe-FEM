# $Id: makefile,v 1.2 2002-09-30 15:26:36 sawimme Exp $
# INSTALLATION INSTRUCTIONS
#
# (0) copy makefile.template to a local makefile.
# (1) set HOME_DIR to the full path to this makefile.
# (2) set MACRO_DIR to the full path to directory
#     containing the machine specific compiler macros.
# (3) select ARCH from those in the macros directory.
# (4) select installed "build options".
# (5) "make init" to (re-)initialize headers and file
#     dependencies.
# (6) "make build" to (re-)compiler binaries.

# full path to current directory
HOME_DIR = $(TAHOE)/contrib/MakeCSE

# working directories:
OBJ_DIR = $(HOME_DIR)/obj
INC_DIR = $(HOME_DIR)/inc
SRC_DIR = $(HOME_DIR)/src
LIB_DIR = $(HOME_DIR)/lib

# build architecture
ARCH = SUN
# see macros directory for choices $(ARCH).macros

# default target name - root of library name or executable name
TARGET  = MakeCSE
LIBRARY = $(HOME_DIR)/lib/lib$(TARGET).a

# path to main, needed if primary target is an executable (not 
# a library)
MAIN    = $(SRC_DIR)/main.o

# target of "make build", i.e., the library $(LIBRARY) or
# the executable $(TARGET)
BUILD   = $(TARGET)

# full path to compiler macros
MACRO_DIR = $(HOME_DIR)/../../macros

############################### build options ###########################
# uncomment macros for installed models

##### toolbox
#TB_DIR      = $(HOME_DIR)/../toolbox
#TB_INC      = -I$(TB_DIR)/inc
#TB_CFLAGS   = $(TB_INC)
#TB_LIB_NAME = toolbox
#TB_LIB      = $(TB_DIR)/lib/lib$(TB_LIB_NAME).a
#TB_LFLAGS   = -L$(TB_DIR)/lib -l$(TB_LIB_NAME)

##### SEACAS/ACCESS software
AX_INC    = -I$(ACCESS)/inc
AX_CFLAGS = $(AX_INC) -D__ACCESS__
AX_LFLAGS = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf
############################### build options ###########################

# include common target definitions
include $(MACRO_DIR)/make.targets
