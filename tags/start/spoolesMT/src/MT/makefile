# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MT.QRfactorMT.o \
	MT.QRsolveMT.o \
	MT.factorMT.o \
	MT.mvm.o \
	MT.solveMT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MT.QRfactorMT.c \
	MT.QRsolveMT.c \
	MT.factorMT.c \
	MT.mvm.c \
	MT.solveMT.c

# headers
HDR = \
	MT.spoolesMT.h
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





