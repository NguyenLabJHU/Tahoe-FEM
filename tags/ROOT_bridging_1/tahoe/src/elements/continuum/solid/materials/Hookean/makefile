# $Id: makefile,v 1.5 2002-06-27 01:22:58 hspark Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FDHookeanMatT.o \
	HookeanMatT.o \
	SSHookeanMatT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FDHookeanMatT.cpp \
	HookeanMatT.cpp \
	SSHookeanMatT.cpp

# headers
HDR = \
	FDHookeanMatT.h \
	HookeanMatT.h \
	SSHookeanMatT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	Cubic \
	1D \
	KStV
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
