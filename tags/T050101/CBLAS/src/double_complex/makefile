# $Id: makefile,v 1.4 2001-02-22 23:40:25 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	dzasum.o \
	dznrm2.o \
	izamax.o \
	zaxpy.o \
	zcopy.o \
	zdotc.o \
	zgemv.o \
	zgerc.o \
	zhemv.o \
	zher2.o \
	zmyblas2.o \
	zscal.o \
	ztrsv.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	dzasum.c \
	dznrm2.c \
	izamax.c \
	zaxpy.c \
	zcopy.c \
	zdotc.c \
	zgemv.c \
	zgerc.c \
	zhemv.c \
	zher2.c \
	zmyblas2.c \
	zscal.c \
	ztrsv.c

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
