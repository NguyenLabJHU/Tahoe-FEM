include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MakeCSE_IOManager.o \
	MakeCSE_FEManager.o \
	MakeCSE_ElementBaseT.o \
	Quad2Tri.o \
	MakeCSE_CSEBaseT.o \
	MakeCSE.o \
	GlobalEdgeFinderT.o \
	NodeManagerPrimitive.o \
	CSEConstants.o \
	InteractiveIOManagerT.o \
	main.o \
	ParameterFileManagerT.o \
	MakeCSE_ExecutionT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GlobalEdgeFinderT.cpp \
	MakeCSE_ElementBaseT.cpp \
	MakeCSE_CSEBaseT.cpp \
	MakeCSE_FEManager.cpp \
	MakeCSE.cpp \
	MakeCSE_IOManager.cpp \
	NodeManagerPrimitive.cpp \
	Quad2Tri.cpp \
	CSEConstants.cpp \
	InteractiveIOManagerT.cpp \
	ParameterFileManagerT.cpp \
	main.cpp \
	MakeCSE_ExecutionT.cpp

# headers
HDR = \
	GlobalEdgeFinderT.h \
	MakeCSE_ElementBaseT.h \
	MakeCSE_CSEBaseT.h \
	MakeCSE_FEManager.h \
	MakeCSE.h \
	MakeCSE_IOManager.h \
	NodeManagerPrimitive.h \
	Quad2Tri.h \
	CSEConstants.h \
	InteractiveIOManagerT.h \
	ParameterFileManagerT.h \
	MakeCSE_ExecutionT.h \
	ClockT.h

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
