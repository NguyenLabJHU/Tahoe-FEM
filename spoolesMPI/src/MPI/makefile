# $Id: makefile,v 1.3 2001-02-22 23:43:28 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MPI.aggListMPI.o \
	MPI.colmapMPI.o \
	MPI.DenseMtx_gather.o \
	MPI.DenseMtx_scatterAdd.o \
	MPI.ETree_Bcast.o \
	MPI.factorMPI.o \
	MPI.fullAdjMPI.o \
	MPI.Graph_Bcast.o \
	MPI.IVallgather.o \
	MPI.IV_Bcast.o \
	MPI.IVLallgather.o \
	MPI.IVL_alltoall.o \
	MPI.IVL_Bcast.o \
	MPI.makeSendRecvIVLs.o \
	MPI.MMM.o \
	MPI.postProcess.o \
	MPI.rowmapMPI.o \
	MPI.solveMPI.o \
	MPI.splitDenseMtx.o \
	MPI.splitFrontMtx.o \
	MPI.splitInpMtx.o \
	MPI.splitPencil.o \
	MPI.symbfacMPI.o \
	MPI.utilities.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MPI.aggListMPI.c \
	MPI.colmapMPI.c \
	MPI.DenseMtx_gather.c \
	MPI.DenseMtx_scatterAdd.c \
	MPI.ETree_Bcast.c \
	MPI.factorMPI.c \
	MPI.fullAdjMPI.c \
	MPI.Graph_Bcast.c \
	MPI.IVallgather.c \
	MPI.IV_Bcast.c \
	MPI.IVLallgather.c \
	MPI.IVL_alltoall.c \
	MPI.IVL_Bcast.c \
	MPI.makeSendRecvIVLs.c \
	MPI.MMM.c \
	MPI.postProcess.c \
	MPI.rowmapMPI.c \
	MPI.solveMPI.c \
	MPI.splitDenseMtx.c \
	MPI.splitFrontMtx.c \
	MPI.splitInpMtx.c \
	MPI.splitPencil.c \
	MPI.symbfacMPI.c \
	MPI.utilities.c

# headers
HDR = \
	MPI.spoolesMPI.h
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
