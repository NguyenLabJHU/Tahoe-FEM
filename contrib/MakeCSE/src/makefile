include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	MakeCSEIOManager.o \
	FEManager.o \
	ElementBaseT.o \
	Quad2Tri.o \
	CSEBaseT.o \
	MakeCSE.o \
	GlobalEdgeFinderT.o \
	NodeManagerPrimitive.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GlobalEdgeFinderT.cpp \
	ElementBaseT.cpp \
	CSEBaseT.cpp \
	FEManager.cpp \
	MakeCSE.cpp \
	MakeCSEIOManager.cpp \
	NodeManagerPrimitive.cpp \
	Quad2Tri.cpp \
	IOManager.cpp \
	main.cpp

# headers
HDR = \
	GlobalEdgeFinderT.h \
	ElementBaseT.h \
	CSEBaseT.h \
	FEManager.h \
	MakeCSE.h \
	MakeCSEIOManager.h \
	NodeManagerPrimitive.h \
	Quad2Tri.h \
	IOManager.h \
	ClockT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	geometry
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
