# $Id: makefile,v 1.5 2005-08-12 05:31:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	abort_.o \
	cabs.o \
	derf_.o \
	derfc_.o \
	erf_.o \
	erfc_.o \
	exit_.o \
	f77vers.o \
	getarg_.o \
	getenv_.o \
	iargc_.o \
	i77vers.o \
	sig_die.o \
	signal_.o \
	s_paus.o \
	s_rnge.o \
	s_stop.o \
	system_.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	abort_.c \
	cabs.c \
	derf_.c \
	derfc_.c \
	erf_.c \
	erfc_.c \
	exit_.c \
	f77vers.c \
	getarg_.c \
	getenv_.c \
	iargc_.c \
	i77vers.c \
	sig_die.c \
	signal_.c \
	s_paus.c \
	s_rnge.c \
	s_stop.c \
	system_.c

# headers
HDR = \
	fp.h \
	f2c.h \
	libf2cCWPrefix.h \
	signal1.h
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
