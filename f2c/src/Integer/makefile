# $Id: makefile,v 1.5 2005-08-12 05:31:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	i_abs.o \
	i_dim.o \
	i_dnnt.o \
	i_indx.o \
	i_len.o \
	i_mod.o \
	i_nint.o \
	i_sign.o \
	lbitbits.o \
	lbitshft.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	i_abs.c \
	i_dim.c \
	i_dnnt.c \
	i_indx.c \
	i_len.c \
	i_mod.c \
	i_nint.c \
	i_sign.c \
	lbitbits.c \
	lbitshft.c

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
