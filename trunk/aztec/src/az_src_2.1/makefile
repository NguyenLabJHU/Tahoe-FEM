# $Id: makefile,v 1.2 2003-02-28 01:19:17 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	az_bilu.o \
	az_cg.o \
	az_cgs.o \
	az_cgstab.o \
	az_check.o \
	az_colperm.o \
	az_comm.o \
	az_converge.o \
	az_dd_overlap.o \
	az_dgemv2.o \
	az_dgemv3.o \
	az_domain_decomp.o \
	az_fix_pt.o \
	az_flop_cnt.o \
	az_fortran_wrap.o \
	az_f_util.o \
	az_gmres.o \
	az_gmresr.o \
	az_gsumd_puma.o \
	az_icc.o \
	az_ilu_util.o \
	az_ilut.o \
	az_interface.o \
	az_lu_y12.o \
	az_matrix_util.o \
	az_matvec_mult.o \
	az_old_matvec_mult.o \
	az_pad_utils.o \
	az_poly.o \
	az_precond.o \
	az_qmrcgs.o \
	az_reorder.o \
	az_rilu.o \
	az_scaling.o \
	az_solve.o \
	az_sort.o \
	az_subdomain_solver.o \
	az_tools.o \
	az_util.o \
	md_timer_generic.o \
	md_wrap.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	az_bilu.c \
	az_cg.c \
	az_cgs.c \
	az_cgstab.c \
	az_check.c \
	az_colperm.c \
	az_comm.c \
	az_converge.c \
	az_dd_overlap.c \
	az_dgemv2.c \
	az_dgemv3.c \
	az_domain_decomp.c \
	az_fix_pt.c \
	az_flop_cnt.c \
	az_fortran_wrap.c \
	az_f_util.c \
	az_gmres.c \
	az_gmresr.c \
	az_gsumd_puma.c \
	az_icc.c \
	az_ilu_util.c \
	az_ilut.c \
	az_interface.c \
	az_lu_y12.c \
	az_matrix_util.c \
	az_matvec_mult.c \
	az_old_matvec_mult.c \
	az_pad_utils.c \
	az_poly.c \
	az_precond.c \
	az_qmrcgs.c \
	az_reorder.c \
	az_rilu.c \
	az_scaling.c \
	az_solve.c \
	az_sort.c \
	az_subdomain_solver.c \
	az_tools.c \
	az_util.c \
	md_timer_generic.c \
	md_wrap.c

# headers
HDR = \
	az_aztec.h \
	az_aztec_defs.h \
	az_aztecf.h \
	az_blas_wrappers.h \
	az_f77func.h \
	az_lapack_wrappers.h \
	az_y12m_wrappers.h

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
