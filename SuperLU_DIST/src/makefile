# $Id: makefile,v 1.1.1.1 2005-07-30 18:04:14 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# headers
HDR = \
Cnames.h \
colamd.h \
dcomplex.h \
machines.h \
superlu_ddefs.h \
superlu_defs.h \
superlu_zdefs.h \
supermatrix.h \
util_dist.h

HDR_LINK = $(HDR:.h=.h_link)

# sources
SRC_F = \
mc21.f \
mc64ad.f

SRC_C = \
GetDiagU.c \
colamd.c \
comm.c \
dcomplex.c \
ddistribute.c \
dgsequ.c \
dlamch.c \
dlangs.c \
dlaqgs.c \
dldperm.c \
dmemory.c \
dmyblas2.c \
dreadhb.c \
dsp_blas2.c \
dsp_blas3.c \
dutil.c \
etree.c \
get_perm_c.c \
lsame.c \
memory.c \
mmd.c \
pddistribute.c \
pdgsequ.c \
pdgsmv.c \
pdgsmv_AXglobal.c \
pdgsrfs.c \
pdgsrfs_ABXglobal.c \
pdgssvx.c \
pdgssvx_ABglobal.c \
pdgstrf_X1.c \
pdgstrf_irecv.c \
pdgstrs.c \
pdgstrs1.c \
pdgstrs_Bglobal.c \
pdgstrs_Bglobal_Bsend.c \
pdgstrs_lsum.c \
pdlangs.c \
pdlaqgs.c \
pdutil.c \
pxerbla.c \
pzdistribute.c \
pzgsequ.c \
pzgsmv.c \
pzgsmv_AXglobal.c \
pzgsrfs.c \
pzgsrfs_ABXglobal.c \
pzgssvx.c \
pzgssvx_ABglobal.c \
pzgstrf_irecv.c \
pzgstrs.c \
pzgstrs1.c \
pzgstrs_Bglobal.c \
pzgstrs_lsum.c \
pzlangs.c \
pzlaqgs.c \
pzutil.c \
slamch.c \
sp_colorder.c \
sp_ienv.c \
superlu_grid.c \
superlu_timer.c \
symbfact.c \
util.c \
xerbla.c \
zdistribute.c \
zdistribute_mark.c \
zgsequ.c \
zlangs.c \
zlaqgs.c \
zldperm.c \
zmemory.c \
zmyblas2.c \
zreadhb.c \
zsp_blas2.c \
zsp_blas3.c \
zutil.c

# objects
OBJ_F = $(SRC_F:.f=.o)
OBJ_C = $(SRC_C:.c=.o)
OBJ = $(OBJ_F) $(OBJ_C)
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

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
