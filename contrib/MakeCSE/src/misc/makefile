# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	StringT.o \
	dMatrixT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	StringT.cpp \
	dMatrixT.cpp

# headers
HDR = \
	StringT.h \
	dMatrixT.h \
	pArrayT.h \
	AutoFill2DT.h \
	nMatrixT.h \
	iosfwd.h \
	ios_fwd_decl.h \
	Constants.h \
	Environment.h \
	ExceptionCodes.h

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
