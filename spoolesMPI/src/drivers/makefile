# $Id: makefile,v 1.7 2004-09-07 06:41:58 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	LU_MPI_driver.o \
	LU_MPI_driver_init.o \
	LU_MPI_driver_factorize.o \
	LU_MPI_driver_solve.o \
	LU_MPI_driver_free.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	LU_MPI_driver.c \
	LU_MPI_driver_init.c \
	LU_MPI_driver_factorize.c \
	LU_MPI_driver_solve.c \
	LU_MPI_driver_free.c

# headers
HDR = LU_MPI_driver.h LU_MPI_driver_int.h
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
