# $Id: makefile,v 1.4 2001-02-22 23:40:05 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	az_cg.o \
	az_cgs.o \
	az_cgstab.o \
	az_check.o \
	az_comm.o \
	az_converge.o \
	az_dgemv2.o \
	az_dgemv3.o \
	az_domain_decomp.o \
	az_domain_setup.o \
	az_flop_cnt.o \
	az_fortran_wrap.o \
	az_gmres.o \
	az_gsumd_puma.o \
	az_icc.o \
	az_ilu_util.o \
	az_ilu_vbr.o \
	az_interface.o \
	az_lu_y12.o \
	az_matrix_util.o \
	az_msr_ilu.o \
	az_poly.o \
	az_precond.o \
	az_qmrcgs.o \
	az_scaling.o \
	az_solve.o \
	az_sort.o \
	az_sparax.o \
	az_symmlq.o \
	az_tools.o \
	az_util.o \
	md_timer_generic.o \
	md_wrap_c.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	az_cg.c \
	az_cgs.c \
	az_cgstab.c \
	az_check.c \
	az_comm.c \
	az_converge.c \
	az_dgemv2.c \
	az_dgemv3.c \
	az_domain_decomp.c \
	az_domain_setup.c \
	az_flop_cnt.c \
	az_fortran_wrap.c \
	az_gmres.c \
	az_gsumd_puma.c \
	az_icc.c \
	az_ilu_util.c \
	az_ilu_vbr.c \
	az_interface.c \
	az_lu_y12.c \
	az_matrix_util.c \
	az_msr_ilu.c \
	az_poly.c \
	az_precond.c \
	az_qmrcgs.c \
	az_scaling.c \
	az_solve.c \
	az_sort.c \
	az_sparax.c \
	az_symmlq.c \
	az_tools.c \
	az_util.c \
	md_timer_generic.c \
	md_wrap_c.c

# headers
HDR = \
	az_aztec_defs.h \
	az_aztec.h
HDR_LINK = $(HDR:.h=.h_link)

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
