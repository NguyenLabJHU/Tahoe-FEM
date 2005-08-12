# $Id: makefile,v 1.5 2005-08-12 05:31:57 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	r_abs.o \
	r_acos.o \
	r_asin.o \
	r_atan.o \
	r_atn2.o \
	r_cnjg.o \
	r_cos.o \
	r_cosh.o \
	r_dim.o \
	r_exp.o \
	r_imag.o \
	r_int.o \
	r_lg10.o \
	r_log.o \
	r_mod.o \
	r_nint.o \
	r_sign.o \
	r_sin.o \
	r_sinh.o \
	r_sqrt.o \
	r_tan.o \
	r_tanh.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	r_abs.c \
	r_acos.c \
	r_asin.c \
	r_atan.c \
	r_atn2.c \
	r_cnjg.c \
	r_cos.c \
	r_cosh.c \
	r_dim.c \
	r_exp.c \
	r_imag.c \
	r_int.c \
	r_lg10.c \
	r_log.c \
	r_mod.c \
	r_nint.c \
	r_sign.c \
	r_sin.c \
	r_sinh.c \
	r_sqrt.c \
	r_tan.c \
	r_tanh.c

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
