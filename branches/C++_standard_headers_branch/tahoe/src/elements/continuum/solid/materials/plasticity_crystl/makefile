include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

##################### (1) objects
#
# (1.1) If there are objects in this directory:
#    (1) define the OBJ macro as the list of object files
#    (2) uncomment the OBJ_LINK and DEP macros here
#    (3) comment the OBJ, OBJ_LINK, and DEP macros in the
#        section (1.1).
#
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
	GradientTools_C.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#
# (1.2) If there are not object in this directory:
#    (1) comment the OBJ, OBJ_LINK, and DEP definitions in section (1.1)
#    (2) uncomment the OBJ, OBJ_LINK, and DEP definitions below
#
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
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
	GradientTools_C.cpp

##################### (3) headers
#
# (3) If there are header files in this directory:
#    (1) list the header files in the HDR macro
#    (2) uncomment the HDR_LINK
# 
# If there are no header files in this directory comment both the
# HDR and HDR_LINK macro   
#
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
	GradientTools_C.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
#SUB_DIR = 
#subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
SUB_DIR =
subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
