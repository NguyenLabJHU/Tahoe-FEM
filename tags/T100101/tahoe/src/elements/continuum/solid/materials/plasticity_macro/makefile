include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	EVPFDBaseT.o \
	HyperEVP3D.o \
	HyperEVP2D.o \
	BCJHypo3D.o  \
	BCJHypo2D.o  \
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
	BCJKineticEqn.cpp \
	SimplePowerLaw.cpp

# headers
HDR = \
	EVPFDBaseT.h \
	HyperEVP3D.h \
	HyperEVP2D.h \
	BCJHypo3D.h  \
	BCJHypo2D.h  \
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
