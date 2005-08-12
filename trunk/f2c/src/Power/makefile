# $Id: makefile,v 1.5 2005-08-12 05:31:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	pow_ci.o \
	pow_dd.o \
	pow_di.o \
	pow_hh.o \
	pow_ii.o \
	pow_ri.o \
	pow_zi.o \
	pow_zz.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	pow_ci.c \
	pow_dd.c \
	pow_di.c \
	pow_hh.c \
	pow_ii.c \
	pow_ri.c \
	pow_zi.c \
	pow_zz.c

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
