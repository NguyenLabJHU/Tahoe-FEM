# $Id: makefile,v 1.4 2001-02-22 23:46:59 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DRSolver.o \
	ExpCD_DRSolver.o \
	iNLSolver_LS.o \
	LinearSolver.o \
	NLK0Solver.o \
	NLSolver.o \
	NLSolver_LS.o \
	NLSolverX.o \
	PCGSolver_LS.o \
	SolverT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DRSolver.cpp \
	ExpCD_DRSolver.cpp \
	iNLSolver_LS.cpp \
	LinearSolver.cpp \
	NLK0Solver.cpp \
	NLSolver.cpp \
	NLSolver_LS.cpp \
	NLSolverX.cpp \
	PCGSolver_LS.cpp \
	SolverT.cpp

# headers
HDR = \
	DRSolver.h \
	ExpCD_DRSolver.h \
	iNLSolver_LS.h \
	LinearSolver.h \
	NLK0Solver.h \
	NLSolver.h \
	NLSolver_LS.h \
	NLSolverX.h \
	PCGSolver_LS.h \
	SolverT.h
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
