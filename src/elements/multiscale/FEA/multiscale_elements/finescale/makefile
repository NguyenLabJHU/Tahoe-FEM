#/* $Id: makefile,v 1.3 2002-12-21 01:42:01 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VMS_BCJT.o \
	VMS_EZT.o	\
	VMS_EZ2T.o	\
	VMS_EZ3T.o	\
	VMS_EZ4T.o	\
	VMS_EZ5T.o	


OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	VMS_BCJT.cpp	\
	VMS_EZT.cpp		\
	VMS_EZ2T.cpp		\
	VMS_EZ3T.cpp		\
	VMS_EZ4T.cpp		\
	VMS_EZ5T.cpp		

# headers
HDR = \
	FineScaleT.h 	\
	VMS_BCJT.h 	\
  VMS_EZT.h	\
  VMS_EZ2T.h	\
  VMS_EZ3T.h \
  VMS_EZ4T.h \
  VMS_EZ5T.h 
	
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
