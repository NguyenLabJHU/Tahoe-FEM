# $Id: makefile,v 1.1.1.1 2005-07-31 01:06:47 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# headers
HDR = \
Cnames.h \
colamd.h \
csp_defs.h \
dcomplex.h \
dsp_defs.h \
scomplex.h \
ssp_defs.h \
supermatrix.h \
util.h \
zsp_defs.h \

HDR_LINK = $(HDR:.h=.h_link)

# sources
SRC = \
ccolumn_bmod.c \
ccolumn_dfs.c \
ccopy_to_ucol.c \
cgscon.c \
cgsequ.c \
cgsrfs.c \
cgssv.c \
cgssvx.c \
cgstrf.c \
cgstrs.c \
clacon.c \
clangs.c \
claqgs.c \
cmemory.c \
cmyblas2.c \
colamd.c \
cpanel_bmod.c \
cpanel_dfs.c \
cpivotL.c \
cpivotgrowth.c \
cpruneL.c \
creadhb.c \
csnode_bmod.c \
csnode_dfs.c \
csp_blas2.c \
csp_blas3.c \
cutil.c \
dGetDiagU.c \
dcolumn_bmod.c \
dcolumn_dfs.c \
dcomplex.c \
dcopy_to_ucol.c \
dgscon.c \
dgsequ.c \
dgsrfs.c \
dgssv.c \
dgssvx.c \
dgstrf.c \
dgstrs.c \
dgstrsL.c \
dlacon.c \
dlamch.c \
dlangs.c \
dlaqgs.c \
dmemory.c \
dmyblas2.c \
dpanel_bmod.c \
dpanel_dfs.c \
dpivotL.c \
dpivotgrowth.c \
dpruneL.c \
dreadhb.c \
dsnode_bmod.c \
dsnode_dfs.c \
dsp_blas2.c \
dsp_blas3.c \
dutil.c \
dzsum1.c \
get_perm_c.c \
heap_relax_snode.c \
icmax1.c \
izmax1.c \
lsame.c \
memory.c \
mmd.c \
relax_snode.c \
scolumn_bmod.c \
scolumn_dfs.c \
scomplex.c \
scopy_to_ucol.c \
scsum1.c \
sgscon.c \
sgsequ.c \
sgsrfs.c \
sgssv.c \
sgssvx.c \
sgstrf.c \
sgstrs.c \
slacon.c \
slamch.c \
slangs.c \
slaqgs.c \
smemory.c \
smyblas2.c \
sp_coletree.c \
sp_ienv.c \
sp_preorder.c \
spanel_bmod.c \
spanel_dfs.c \
spivotL.c \
spivotgrowth.c \
spruneL.c \
sreadhb.c \
ssnode_bmod.c \
ssnode_dfs.c \
ssp_blas2.c \
ssp_blas3.c \
superlu_timer.c \
sutil.c \
util.c \
xerbla.c \
zcolumn_bmod.c \
zcolumn_dfs.c \
zcopy_to_ucol.c \
zgscon.c \
zgsequ.c \
zgsrfs.c \
zgssv.c \
zgssvx.c \
zgstrf.c \
zgstrs.c \
zlacon.c \
zlangs.c \
zlaqgs.c \
zmemory.c \
zmyblas2.c \
zpanel_bmod.c \
zpanel_dfs.c \
zpivotL.c \
zpivotgrowth.c \
zpruneL.c \
zreadhb.c \
zsnode_bmod.c \
zsnode_dfs.c \
zsp_blas2.c \
zsp_blas3.c \
zutil.c

# objects
OBJ = $(SRC:.c=.o)
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
