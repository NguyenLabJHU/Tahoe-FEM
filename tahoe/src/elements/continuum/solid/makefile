# $Id: makefile,v 1.21 2004-06-26 18:39:04 paklein Exp $
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
	SolidElementT.o \
	FiniteStrainT.o \
	LocalizerT.o \
	SimoFiniteStrainT.o \
	SimoQ1P0.o \
	SmallStrainT.o \
	TotalLagrangianT.o \
	UpdatedLagrangianT.o \
	UpLagr_ExternalFieldT.o \
	UpLagAdaptiveT.o \
	SmallStrainAxiT.o \
	FiniteStrainAxiT.o \
	TotalLagrangianAxiT.o \
	UpdatedLagrangianAxiT.o \
	SimoQ1P0Axi.o \
	SimoQ1P0_inv.o \
	SimoQ1P0Axi_inv.o

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
	SolidElementT.cpp \
	FiniteStrainT.cpp \
	LocalizerT.cpp \
	SimoFiniteStrainT.cpp \
	SimoQ1P0.cpp \
	SmallStrainT.cpp \
	TotalLagrangianT.cpp \
	UpdatedLagrangianT.cpp \
	UpLagr_ExternalFieldT.cpp \
	UpLagAdaptiveT.cpp \
	SmallStrainAxiT.cpp \
	FiniteStrainAxiT.cpp \
	TotalLagrangianAxiT.cpp \
	UpdatedLagrangianAxiT.cpp \
	SimoQ1P0Axi.cpp \
	SimoQ1P0_inv.cpp \
	SimoQ1P0Axi_inv.cpp

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
	SolidElementT.h \
	FiniteStrainT.h \
	LocalizerT.h \
	SimoFiniteStrainT.h \
	SimoQ1P0.h \
	SmallStrainT.h \
	TotalLagrangianT.h \
	UpdatedLagrangianT.h \
	UpLagr_ExternalFieldT.h \
	UpLagAdaptiveT.h \
	SmallStrainAxiT.h \
	FiniteStrainAxiT.h \
	TotalLagrangianAxiT.h \
	UpdatedLagrangianAxiT.h \
	SimoQ1P0Axi.h \
	SimoQ1P0_inv.h \
	SimoQ1P0Axi_inv.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
SUB_DIR = \
	materials \
	meshfree 

subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
#SUB_DIR =
#subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
