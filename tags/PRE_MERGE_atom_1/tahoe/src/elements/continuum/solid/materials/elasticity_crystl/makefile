include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CrystalElast.o \
	FDCrystalElast.o \
	CrystalElastLat.o \
	CrystalElastMat.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
        CrystalElast.cpp \
        FDCrystalElast.cpp \
        CrystalElastLat.cpp \
        CrystalElastMat.cpp 

# headers
HDR = \
        CrystalElast.h \
        FDCrystalElast.h \
        CrystalElastLat.h \
        CrystalElastMat.h 

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
