# $Id: makefile,v 1.5 2003-02-28 01:19:16 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	la_dgeco.o \
	la_dgedi.o \
	la_dgefa.o \
	la_dgetf2.o \
	la_dgetrf.o \
	la_dgetri.o \
	la_dgetrs.o \
	la_dlaic1.o \
	la_dlamch.o \
	la_dlaswp.o \
	la_dpotf2.o \
	la_dpotrf.o \
	la_dtrti2.o \
	la_dtrtri.o \
	la_ilaenv.o \
	la_slamch.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	la_dgeco.c \
	la_dgedi.c \
	la_dgefa.c \
	la_dgetf2.c \
	la_dgetrf.c \
	la_dgetri.c \
	la_dgetrs.c \
	la_dlaic1.c \
	la_dlamch.c \
	la_dlaswp.c \
	la_dpotf2.c \
	la_dpotrf.c \
	la_dtrti2.c \
	la_dtrtri.c \
	la_ilaenv.c \
	la_slamch.c

# headers
HDR = \
	linpack.h
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
