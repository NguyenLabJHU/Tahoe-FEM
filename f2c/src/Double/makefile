# $Id: makefile,v 1.5 2005-08-12 05:31:54 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	d_abs.o \
	d_acos.o \
	d_asin.o \
	d_atan.o \
	d_atn2.o \
	d_cnjg.o \
	d_cos.o \
	d_cosh.o \
	d_dim.o \
	d_exp.o \
	d_imag.o \
	d_int.o \
	d_lg10.o \
	d_log.o \
	d_mod.o \
	d_nint.o \
	d_prod.o \
	d_sign.o \
	d_sin.o \
	d_sinh.o \
	d_sqrt.o \
	d_tan.o \
	d_tanh.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	d_abs.c \
	d_acos.c \
	d_asin.c \
	d_atan.c \
	d_atn2.c \
	d_cnjg.c \
	d_cos.c \
	d_cosh.c \
	d_dim.c \
	d_exp.c \
	d_imag.c \
	d_int.c \
	d_lg10.c \
	d_log.c \
	d_mod.c \
	d_nint.c \
	d_prod.c \
	d_sign.c \
	d_sin.c \
	d_sinh.c \
	d_sqrt.c \
	d_tan.c \
	d_tanh.c

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
