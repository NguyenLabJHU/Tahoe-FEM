# $Id: makefile,v 1.1 2002-08-01 23:34:56 rdorgan Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	LocalJ2SSNonlinHard.o\
	LocalJ2SSNonlinHard2D.o\
	GradJ2SSNonlinHard.o\
	GradJ2SSNonlinHard2D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	LocalJ2SSNonlinHard.cpp\
	LocalJ2SSNonlinHard2D.cpp\
	GradJ2SSNonlinHard.cpp\
	GradJ2SSNonlinHard2D.cpp

# headers
HDR = \
	LocalJ2SSNonlinHard.h\
	LocalJ2SSNonlinHard2D.h\
	GradJ2SSNonlinHard.h\
	GradJ2SSNonlinHard2D.h

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
