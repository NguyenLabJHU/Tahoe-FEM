# $Id: makefile,v 1.4 2001-02-22 23:40:04 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	y12cck.o \
	y12mae.o \
	y12maf.o \
	y12mbe.o \
	y12mbf.o \
	y12mce.o \
	y12mcf.o \
	y12mde.o \
	y12mdf.o \
	y12mfe.o \
	y12mge.o \
	y12mhe.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	y12cck.c \
	y12mae.c \
	y12maf.c \
	y12mbe.c \
	y12mbf.c \
	y12mce.c \
	y12mcf.c \
	y12mde.c \
	y12mdf.c \
	y12mfe.c \
	y12mge.c \
	y12mhe.c

# headers
HDR = \
	y12m.h
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
