# $Id: makefile,v 1.4 2001-02-22 23:40:54 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	h_abs.o \
	h_dim.o \
	h_dnnt.o \
	h_indx.o \
	h_len.o \
	h_mod.o \
	h_nint.o \
	h_sign.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	h_abs.c \
	h_dim.c \
	h_dnnt.c \
	h_indx.c \
	h_len.c \
	h_mod.c \
	h_nint.c \
	h_sign.c

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
DEPEND = /dev/null
include $(DEPEND)
