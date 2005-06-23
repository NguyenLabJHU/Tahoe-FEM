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
	EVPFDBaseT.o \
	HyperEVP3D.o \
	HyperEVP2D.o \
	BCJHypo3D.o  \
	BCJHypo2D.o  \
	BCJHypoIsoDamageKE3D.o \
	BCJHypoIsoDamageYC3D.o \
	BCJHypoIsoDamageKE2D.o \
	BCJHypoIsoDamageYC2D.o \
	VoidGrowthModelImp.o \
	BCJKineticEqn.o \
	SimplePowerLaw.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#
# (1.2) If there are not object in this directory:
#    (1) comment the OBJ, OBJ_LINK, and DEP definitions in section (1.1)
#    (2) uncomment the OBJ, OBJ_LINK, and DEP definitions below
#
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
SRC = \
	EVPFDBaseT.cpp \
	HyperEVP3D.cpp \
	HyperEVP2D.cpp \
	BCJHypo3D.cpp  \
	BCJHypo2D.cpp  \
	BCJHypoIsoDamageKE3D.cpp \
	BCJHypoIsoDamageYC3D.cpp \
	BCJHypoIsoDamageKE2D.cpp \
	BCJHypoIsoDamageYC2D.cpp \
	VoidGrowthModelImp.cpp \
	BCJKineticEqn.cpp \
	SimplePowerLaw.cpp

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
	EVPFDBaseT.h \
	HyperEVP3D.h \
	HyperEVP2D.h \
	BCJHypo3D.h  \
	BCJHypo2D.h  \
	BCJHypoIsoDamageKE3D.h \
	BCJHypoIsoDamageYC3D.h \
	BCJHypoIsoDamageKE2D.h \
	BCJHypoIsoDamageYC2D.h \
	VoidGrowthModelImp.h \
	KineticEqnBase.h \
	BCJKineticEqn.h \
	SimplePowerLaw.h

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
DEPEND = /dev/null
include $(DEPEND)
