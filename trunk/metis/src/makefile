# $Id: makefile,v 1.1.1.1 2001-07-19 06:24:00 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	balance.o \
	bucketsort.o \
	ccgraph.o \
	coarsen.o \
	compress.o \
	debug.o \
	estmem.o \
	fm.o \
	fortran.o \
	frename.o \
	graph.o \
	initpart.o \
	kmetis.o \
	kvmetis.o \
	kwayfm.o \
	kwayrefine.o \
	kwayvolfm.o \
	kwayvolrefine.o \
	match.o \
	mbalance.o \
	mbalance2.o \
	mcoarsen.o \
	memory.o \
	mesh.o \
	meshpart.o \
	mfm.o \
	mfm2.o \
	mincover.o \
	minitpart.o \
	minitpart2.o \
	mkmetis.o \
	mkwayfmh.o \
	mkwayrefine.o \
	mmatch.o \
	mmd.o \
	mpmetis.o \
	mrefine.o \
	mrefine2.o \
	mutil.o \
	myqsort.o \
	ometis.o \
	parmetis.o \
	pmetis.o \
	pqueue.o \
	refine.o \
	separator.o \
	sfm.o \
	srefine.o \
	stat.o \
	subdomains.o \
	timing.o \
	util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	balance.c \
	bucketsort.c \
	ccgraph.c \
	coarsen.c \
	compress.c \
	debug.c \
	estmem.c \
	fm.c \
	fortran.c \
	frename.c \
	graph.c \
	initpart.c \
	kmetis.c \
	kvmetis.c \
	kwayfm.c \
	kwayrefine.c \
	kwayvolfm.c \
	kwayvolrefine.c \
	match.c \
	mbalance.c \
	mbalance2.c \
	mcoarsen.c \
	memory.c \
	mesh.c \
	meshpart.c \
	mfm.c \
	mfm2.c \
	mincover.c \
	minitpart.c \
	minitpart2.c \
	mkmetis.c \
	mkwayfmh.c \
	mkwayrefine.c \
	mmatch.c \
	mmd.c \
	mpmetis.c \
	mrefine.c \
	mrefine2.c \
	mutil.c \
	myqsort.c \
	ometis.c \
	parmetis.c \
	pmetis.c \
	pqueue.c \
	refine.c \
	separator.c \
	sfm.c \
	srefine.c \
	stat.c \
	subdomains.c \
	timing.c \
	util.c

# headers
HDR = \
	defs.h \
	macros.h \
	metis.h \
	proto.h \
	rename.h \
	struct.h
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
