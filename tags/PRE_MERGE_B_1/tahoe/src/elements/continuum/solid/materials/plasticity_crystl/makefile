include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	PolyCrystalMatT.o \
	LocalCrystalPlast.o \
	LocalCrystalPlast2D.o \
	LocalCrystalPlast_C.o \
	LocalCrystalPlastFp.o \
	LocalCrystalPlastFp2D.o \
	LocalCrystalPlastFp_C.o \
	GradCrystalPlast.o \
	GradCrystalPlast2D.o \
	GradCrystalPlastFp.o \
	GradCrystalPlastFp2D.o \
	SlipGeometry.o \
	LatticeOrient.o \
	CrystalElasticity.o \
	SlipKinetics.o \
	PowerLawIKinetics.o \
	PowerLawIIKinetics.o \
	HaasenKinetics.o \
	SlipHardening.o \
	VoceHardening.o \
	HaasenHardening.o \
	VoceGradHardening.o \
	GradientTools.o \
	Utils.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	PolyCrystalMatT.cpp \
	LocalCrystalPlast.cpp \
	LocalCrystalPlast2D.cpp \
	LocalCrystalPlast_C.cpp \
	LocalCrystalPlastFp.cpp \
	LocalCrystalPlastFp2D.cpp \
	LocalCrystalPlastFp_C.cpp \
	GradCrystalPlast.cpp \
	GradCrystalPlast2D.cpp \
	GradCrystalPlastFp.cpp \
	GradCrystalPlastFp2D.cpp \
	SlipGeometry.cpp \
	LatticeOrient.cpp \
	CrystalElasticity.cpp \
	SlipKinetics.cpp \
	PowerLawIKinetics.cpp \
	PowerLawIIKinetics.cpp \
	HaasenKinetics.cpp \
	SlipHardening.cpp \
	VoceHardening.cpp \
	HaasenHardening.cpp \
	VoceGradHardening.cpp \
	GradientTools.cpp \
	Utils.cpp

# headers
HDR = \
	PolyCrystalMatT.h \
	LocalCrystalPlast.h \
	LocalCrystalPlast2D.h \
	LocalCrystalPlast_C.h \
	LocalCrystalPlastFp.h \
	LocalCrystalPlastFp2D.h \
	LocalCrystalPlastFp_C.h \
	GradCrystalPlast.h \
	GradCrystalPlast2D.h \
	GradCrystalPlastFp.h \
	GradCrystalPlastFp2D.h \
	SlipGeometry.h \
	LatticeOrient.h \
	CrystalElasticity.h \
	SlipKinetics.h \
	PowerLawIKinetics.h \
	PowerLawIIKinetics.h \
	HaasenKinetics.h \
	SlipHardening.h \
	VoceHardening.h \
	HaasenHardening.h \
	VoceGradHardening.h \
	GradientTools.h \
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
