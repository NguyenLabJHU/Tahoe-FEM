# $Id: makefile,v 1.2 2003-05-04 22:49:50 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	basement.o \
	bedroom.o \
	closet.o \
	crawl_space.o \
	driveway.o \
	garage.o \
	house.o \
	lawn.o \
	main.o \
	roof.o \
	room.o \
	storm_shelter.o \
	window.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	basement.cpp \
	bedroom.cpp \
	closet.cpp \
	crawl_space.cpp \
	driveway.cpp \
	garage.cpp \
	house.cpp \
	lawn.cpp \
	main.cpp \
	roof.cpp \
	room.cpp \
	storm_shelter.cpp \
	window.cpp

# headers
HDR = \
	basement.h \
	bedroom.h \
	closet.h \
	crawl_space.h \
	driveway.h \
	garage.h \
	house.h \
	lawn.h \
	roof.h \
	room.h \
	storm_shelter.h \
	window.h

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
