# $Id: makefile,v 1.2 2005-04-05 16:04:19 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	LU_MT_driver.o \
	LU_MT_driver_factorize.o \
	LU_MT_driver_free.o \
	LU_MT_driver_init.o \
	LU_MT_driver_solve.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	LU_MT_driver.c \
	LU_MT_driver_factorize.c \
	LU_MT_driver_free.c \
	LU_MT_driver_init.c \
	LU_MT_driver_solve.c

# headers
HDR = 
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








