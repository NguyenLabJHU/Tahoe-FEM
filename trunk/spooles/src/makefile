# $Id: makefile,v 1.4 2005-08-12 05:27:09 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
#OBJ = 
#OBJ_LINK = $(OBJ)
#DEP = $(OBJ:.o=.d)
OBJ = dummy
OBJ_LINK =
DEP = $(MACRO_DIR)/empty

# sources
#SRC = 

# headers
HDR = \
	cfiles.h \
	Make.inc \
	SPOOLES.h \
	timings.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	A2 \
	BKL \
	BPG \
	Chv \
	ChvList \
	ChvManager \
	Coords \
	DenseMtx \
	Drand \
	drivers \
	DSTree \
	DV \
	EGraph \
	ETree \
	FrontMtx \
	GPart \
	Graph \
	Ideq \
	IIheap \
	ILUMtx \
	InpMtx \
	Iter \
	IV \
	IVL \
	I2Ohash \
	Lock \
	misc \
	MSMD \
	Network \
	PatchAndGoInfo \
	Pencil \
	Perm \
	SemiImplMtx \
	SolveMap \
	SubMtx \
	SubMtxList \
	SubMtxManager \
	SymbFac \
	Tree \
	Utilities \
	ZV
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
