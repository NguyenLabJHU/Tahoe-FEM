# $Id: makefile,v 1.5 2005-08-12 05:31:55 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	backspace.o \
	close.o \
	dfe.o \
	dolio.o \
	due.o \
	endfile.o \
	err.o \
	fmt.o \
	fmtlib.o \
	ftell_.o \
	iio.o \
	ilnw.o \
	inquire.o \
	lread.o \
	lwrite.o \
	open.o \
	rdfmt.o \
	rewind.o \
	rsfe.o \
	rsli.o \
	rsne.o \
	sfe.o \
	sue.o \
	typesize.o \
	uio.o \
	util.o \
	wref.o \
	wrtfmt.o \
	wsfe.o \
	wsle.o \
	wsne.o \
	xwsne.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	backspace.c \
	close.c \
	dfe.c \
	dolio.c \
	due.c \
	endfile.c \
	err.c \
	fmt.c \
	fmtlib.c \
	ftell_.c \
	iio.c \
	ilnw.c \
	inquire.c \
	lread.c \
	lwrite.c \
	open.c \
	rdfmt.c \
	rewind.c \
	rsfe.c \
	rsli.c \
	rsne.c \
	sfe.c \
	sue.c \
	typesize.c \
	uio.c \
	util.c \
	wref.c \
	wrtfmt.c \
	wsfe.c \
	wsle.c \
	wsne.c \
	xwsne.c

# headers
HDR = \
	fio.h \
	fmt.h \
	lio.h \
	rawio.h
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
