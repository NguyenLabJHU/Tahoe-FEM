include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
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
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
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

# headers
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
