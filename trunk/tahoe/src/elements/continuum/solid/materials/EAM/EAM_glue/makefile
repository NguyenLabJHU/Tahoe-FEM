# $Id: makefile,v 1.4 2001-02-22 23:46:25 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ErcolessiAdamsAl.o \
	FBD_EAMGlue.o \
	VoterChenAl.o \
	VoterChenCu.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ErcolessiAdamsAl.cpp \
	FBD_EAMGlue.cpp \
	VoterChenAl.cpp \
	VoterChenCu.cpp

# headers
HDR = \
	ErcolessiAdamsAl.h \
	FBD_EAMGlue.h \
	VoterChenAl.h \
	VoterChenCu.h
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
