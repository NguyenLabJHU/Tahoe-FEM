#
# This makefile is divided into four sections:
#    (1) objects
#    (2) sources
#    (3) headers
#    (4) subdirectories
#
# Instructions for how to define these is given in each section below
#
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
	FSDEMatSupportQ1P02DT.o \
	FSDEMatQ1P02DT.o \
	Mech_PK2_Q1P02D.o \
	Mech_TanMod_Q1P02D.o \
	ME_MixedModulus_Q1P02D.o \
	ME_MixedModulus_Q1P02Dspatial.o \
	ME_PK2_Q1P02D.o \
	ME_TanMod_Q1P02D.o \
	Elec_PK2_Q1P02D.o

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
	FSDEMatSupportQ1P02DT.cpp \
	FSDEMatQ1P02DT.cpp \
	Mech_PK2_Q1P02D.cpp \
	Mech_TanMod_Q1P02D.cpp \
	ME_MixedModulus_Q1P02D.cpp \
	ME_MixedModulus_Q1P02Dspatial.cpp \
	ME_PK2_Q1P02D.cpp \
	ME_TanMod_Q1P02D.cpp \
	Elec_PK2_Q1P02D.cpp

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
	FSDEMatSupportQ1P02DT.h \
	FSDEMatSupportQ1P02DT.i.h \
	FSDEMatQ1P02DT.h \
	FSDEMatQ1P02DT.i.h \
	FSDE_common_definesQ1P02D.h \
	FSDE_incQ1P02D.h 

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
