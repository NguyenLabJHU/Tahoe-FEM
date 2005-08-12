# $Id: makefile,v 1.5 2005-08-12 05:31:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	hl_ge.o \
	hl_gt.o \
	hl_le.o \
	hl_lt.o \
	l_ge.o \
	l_gt.o \
	l_le.o \
	l_lt.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	hl_ge.c \
	hl_gt.c \
	hl_le.c \
	hl_lt.c \
	l_ge.c \
	l_gt.c \
	l_le.c \
	l_lt.c

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
