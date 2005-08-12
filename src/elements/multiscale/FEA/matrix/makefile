#/* $Id: makefile,v 1.2 2005-08-12 05:38:48 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FEA_dMatrixT.o 	\
	FEA_dVectorT.o 	\
	FEA_dScalarT.o 	\
	FEA_EquateT.o 	\
	FEA_StackT.o 		\
	FEA_ArraysT.o 	\
	FEA_Access.o 		

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	FEA_dMatrixT.cpp 	\
	FEA_dVectorT.cpp 	\
	FEA_dScalarT.cpp 	\
	FEA_EquateT.cpp 	\
	FEA_StackT.cpp 		\
	FEA_ArraysT.cpp 	\
	FEA_Access.cpp 		

# headers
HDR = \
	FEA_dMatrixT.h 	\
	FEA_dVectorT.h 	\
	FEA_dScalarT.h 	\
	FEA_EquateT.h 	\
	FEA_StackT.h 		\
	FEA_ArraysT.h 	
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
