include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	NLCSolver.o \
	NLCSolver_LS.o \
	NLCSolverWrapper.o \
	NLCSolverWrapperPtr.o \
	NewtonMethodBase.o \
	NewtonMethod.o \
	Utils.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	NLCSolver.cpp \
	NLCSolver_LS.cpp \
	NLCSolverWrapper.cpp \
	NLCSolverWrapperPtr.cpp \
	NewtonMethodBase.cpp \
	NewtonMethod.cpp \
	Utils.cpp

# headers
HDR = \
	NLCSolver.h \
	NLCSolver_LS.h \
	NLCSolverWrapper.h \
	NLCSolverWrapperPtr.h \
	NewtonMethodBase.h \
	NewtonMethod.h \
	Utils.h

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
