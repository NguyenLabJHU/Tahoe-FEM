# $Id: makefile,v 1.5 2005-08-12 05:31:54 paklein Exp $
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
#HDR = 
#HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	Char \
	Complex \
	custom_Misc \
	Double \
	Double_Complex \
	ELF \
	Half_Integer \
	Integer \
	IO \
	Logical \
	Misc \
	Power \
	Quad_Integer \
	Real \
	Time
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
