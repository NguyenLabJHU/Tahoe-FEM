# $Id: makefile,v 1.4 2001-02-22 23:40:25 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	isamax.o \
	sasum.o \
	saxpy.o \
	scopy.o \
	sdot.o \
	sgemv.o \
	sger.o \
	smyblas2.o \
	snrm2.o \
	srot.o \
	sscal.o \
	ssymv.o \
	ssyr2.o \
	strsv.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	isamax.c \
	sasum.c \
	saxpy.c \
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
	strsv.c

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
