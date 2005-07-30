# $Id: makefile,v 1.1.1.1 2005-07-30 16:59:01 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# headers
HDR = \
Cnames.h \
f2c.h \
scomplex.h \
superlu_f2c.h

HDR_LINK = $(HDR:.h=.h_link)

# sources
SRC = \
caxpy.c \
ccopy.c \
cdotc.c \
cgemv.c \
cgerc.c \
cgeru.c \
chemv.c \
cher2.c \
cmyblas2.c \
cscal.c \
ctrsv.c \
dasum.c \
daxpy.c \
dcabs1.c \
dcopy.c \
ddot.c \
dgemm.c \
dgemv.c \
dger.c \
dnrm2.c \
drot.c \
dscal.c \
dsymv.c \
dsyr2.c \
dtrsm.c \
dtrsv.c \
dzasum.c \
dznrm2.c \
icamax.c \
idamax.c \
isamax.c \
izamax.c \
sasum.c \
saxpy.c \
scasum.c \
scnrm2.c \
scopy.c \
sdot.c \
sgemv.c \
sger.c \
smyblas2.c \
snrm2.c \
srot.c \
sscal.c \
ssymv.c \
ssyr2.c \
strsv.c \
zaxpy.c \
zcopy.c \
zdotc.c \
zgemm.c \
zgemv.c \
zgerc.c \
zgeru.c \
zhemv.c \
zher2.c \
zscal.c \
ztrsm.c \
ztrsv.c

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
