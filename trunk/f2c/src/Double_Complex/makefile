# $Id: makefile,v 1.5 2005-08-12 05:31:55 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	z_abs.o \
	z_cos.o \
	z_div.o \
	z_exp.o \
	z_log.o \
	z_sin.o \
	z_sqrt.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	z_abs.c \
	z_cos.c \
	z_div.c \
	z_exp.c \
	z_log.c \
	z_sin.c \
	z_sqrt.c

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
