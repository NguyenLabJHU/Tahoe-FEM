include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GlobalEdgeFinderT.o \
	ElementBaseT.o \
	CSEBaseT.o \
	FEManager.o \
	MakeCSE.o \
	MakeCSEIOManager.o \
	NodeManagerPrimitive.o \
	Quad2Tri.o \
	main.o

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
	ClockT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	array \
	misc \
	ifstream_x \
	geometry \
	dataio
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
