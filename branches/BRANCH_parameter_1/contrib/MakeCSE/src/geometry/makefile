# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GeometryBaseT.o \
	HexahedronT.o \
	LineT.o \
	PentahedronT.o \
	QuadT.o \
	TetrahedronT.o \
	TriT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GeometryBaseT.cpp \
	HexahedronT.cpp \
	LineT.cpp \
	PentahedronT.cpp \
	QuadT.cpp \
	TetrahedronT.cpp \
	TriT.cpp

# headers
HDR = \
	GeometryBaseT.h \
	HexahedronT.h \
	LineT.h \
	PentahedronT.h \
	QuadT.h \
	TetrahedronT.h \
	TriT.h

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
