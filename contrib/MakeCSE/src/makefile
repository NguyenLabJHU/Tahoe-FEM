include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	main.o \
	MakeCSEIOManager.o \
	MakeCSEFEManager.o \
	MakeCSE_ElementBaseT.o \
	Quad2Tri.o \
	MakeCSE_CSEBaseT.o \
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
	MakeCSE_ElementBaseT.cpp \
	MakeCSE_CSEBaseT.cpp \
	MakeCSEFEManager.cpp \
	MakeCSE.cpp \
	MakeCSEIOManager.cpp \
	NodeManagerPrimitive.cpp \
	Quad2Tri.cpp \
	main.cpp

# headers
HDR = \
	GlobalEdgeFinderT.h \
	MakeCSE_ElementBaseT.h \
	MakeCSE_CSEBaseT.h \
	MakeCSEFEManager.h \
	MakeCSE.h \
	MakeCSEIOManager.h \
	NodeManagerPrimitive.h \
	Quad2Tri.h \
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
