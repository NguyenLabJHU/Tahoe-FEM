include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FEDEManagerT.o \
	DEManagerT.o \
	assembly.o \
	GhostParticleT.o \
	particle.o \
	cylinder.o \
	rectangle.o \
	ran.o \
	root6.o \
	zrhqr.o \
	hqr.o \
	balanc.o \
	nrutil.o \
	vec.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	FEDEManagerT.cpp \
	DEManagerT.cpp \
	assembly.cpp \
	GhostParticleT.cpp \
	particle.cpp \
	cylinder.cpp \
	rectangle.cpp \
	ran.cpp \
	root6.cpp \
	zrhqr.cpp \
	hqr.cpp \
	balanc.cpp \
	nrutil.cpp \
	vec.cpp

# headers
HDR = \
	FEDEManagerT.h \
	DEManagerT.h \
	assembly.h \
	boundary.h \
	boundarytgt.h \
	const.h \
	contact.h \
	gradation.h \
	nr.h \
	parameter.h \
	shape.h \
	GhostParticleT.h \
	particle.h \
	cylinder.h \
	rectangle.h \
	ran.h \
	root6.h \
	nrutil.h \
	vec.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
