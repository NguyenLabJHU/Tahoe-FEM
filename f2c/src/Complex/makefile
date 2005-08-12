# $Id: makefile,v 1.5 2005-08-12 05:31:54 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	c_abs.o \
	c_cos.o \
	c_div.o \
	c_exp.o \
	c_log.o \
	c_sin.o \
	c_sqrt.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	c_abs.c \
	c_cos.c \
	c_div.c \
	c_exp.c \
	c_log.c \
	c_sin.c \
	c_sqrt.c

# headers
#HDR = 
#HDR_LINK = $(HDR:.h=.h_link)

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
