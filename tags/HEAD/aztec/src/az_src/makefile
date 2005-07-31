# $Id: makefile,v 1.1.1.1 2001-01-30 20:59:12 paklein Exp $
# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .C .f .F .o .d .h .h_link

# objects
OBJ = \
	az_cg.o \
	az_cgs.o \
	az_cgstab.o \
	az_check.o \
	az_comm.o \
	az_converge.o \
	az_dgemv2.o \
	az_dgemv3.o \
	az_domain_decomp.o \
	az_domain_setup.o \
	az_flop_cnt.o \
	az_fortran_wrap.o \
	az_gmres.o \
	az_gsumd_puma.o \
	az_icc.o \
	az_ilu_util.o \
	az_ilu_vbr.o \
	az_interface.o \
	az_lu_y12.o \
	az_matrix_util.o \
	az_msr_ilu.o \
	az_poly.o \
	az_precond.o \
	az_qmrcgs.o \
	az_scaling.o \
	az_solve.o \
	az_sort.o \
	az_sparax.o \
	az_symmlq.o \
	az_tools.o \
	az_util.o \
	md_timer_generic.o \
	md_wrap_c.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	az_cg.c \
	az_cgs.c \
	az_cgstab.c \
	az_check.c \
	az_comm.c \
	az_converge.c \
	az_dgemv2.c \
	az_dgemv3.c \
	az_domain_decomp.c \
	az_domain_setup.c \
	az_flop_cnt.c \
	az_fortran_wrap.c \
	az_gmres.c \
	az_gsumd_puma.c \
	az_icc.c \
	az_ilu_util.c \
	az_ilu_vbr.c \
	az_interface.c \
	az_lu_y12.c \
	az_matrix_util.c \
	az_msr_ilu.c \
	az_poly.c \
	az_precond.c \
	az_qmrcgs.c \
	az_scaling.c \
	az_solve.c \
	az_sort.c \
	az_sparax.c \
	az_symmlq.c \
	az_tools.c \
	az_util.c \
	md_timer_generic.c \
	md_wrap_c.c

# headers
HDR = \
	az_aztec_defs.h \
	az_aztec.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
#SUB_DIR = 
SUB_DIR = no_subdirs

# process sub-directories as targets:
#subdir_driver: $(SUB_DIR)
subdir_driver:

MACROS = \
	"MAKE = $(MAKE)" \
	"MAKE_OPTS = $(MAKE_OPTS)" \
	"INC_DIR = $(INC_DIR)" \
	"OBJ_DIR = $(OBJ_DIR)" \
	"LIB_DIR = $(LIB_DIR)" \
	"TR = $(TR)" \
	"COMP_CC = $(COMP_CC)" \
	"COMP_C = $(COMP_C)" \
	"COMP_F = $(COMP_F)" \
	"CFLAGS_C = $(CFLAGS_C)" \
	"CFLAGS_CC = $(CFLAGS_CC)" \
	"CFLAGS_F = $(CFLAGS_F)" \
	"LIBRARY = $(LIBRARY)" \
	"AR = $(AR)" \
	"ARFLAGS = $(ARFLAGS)" \
	"RM = $(RM)" \
	"RM_FILES = $(RM_FILES)" \
	"LN = $(LN)" \
	"MAKE_DEPEND = $(MAKE_DEPEND)" \
	"DEPEND_PATH = $(DEPEND_PATH)" \
	"DEPEND = $(DEPEND)"

# debugging
echo_all:
	@ echo "MAKE = $(MAKE)"
	@ echo "MAKE_OPTS = $(MAKE_OPTS)"
	@ echo "INC_DIR = $(INC_DIR)"
	@ echo "OBJ_DIR = $(OBJ_DIR)"
	@ echo "LIB_DIR = $(LIB_DIR)"
	@ echo "TR = $(TR)"
	@ echo "COMP_C = $(COMP_C)"
	@ echo "COMP_CC = $(COMP_CC)"
	@ echo "COMP_F = $(COMP_F)"
	@ echo "CFLAGS_C = $(CFLAGS_C)"
	@ echo "CFLAGS_CC = $(CFLAGS_CC)"
	@ echo "CFLAGS_F = $(CFLAGS_F)"
	@ echo "LIBRARY = $(LIBRARY)"
	@ echo "AR = $(AR)"
	@ echo "ARFLAGS = $(ARFLAGS)"
	@ echo "RM = $(RM)"
	@ echo "RM_FILES = $(RM_FILES)"
	@ echo "LN = $(LN)"
	@ echo "MAKE_DEPEND = $(MAKE_DEPEND)"
	@ echo "DEPEND_PATH = $(DEPEND_PATH)"
	@ echo "DEPEND = $(DEPEND)"

clean:
	-@ $(RM) $(RM_FILES)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = clean"

depend_init: all.depend
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = depend_init"

all.depend: $(DEP)
	@ echo "    cat: $(DEP) > all.depend"; cat $(DEP) > all.depend

headers: $(HDR_LINK)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = headers"

objects: all.depend
	@ $(MAKE) $(MACROS) update_objects
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = objects"

update_objects: $(OBJ)

library:
	@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ_LINK)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = library"

$(SUB_DIR): dummy
	@ echo " cd $@"; \
	cd $@; $(MAKE) $(COMMAND) $(MACROS) "CURR_DIR = $(CURR_DIR)/$@"

dummy:

# suffix rules

.c.o:
	@ echo "    $(COMP_C): "$<; \
	$(COMP_C) $(CFLAGS_C) $<
	@echo "#touch" > $(LIB_DIR)/lib.junk

.cpp.o:
	@ echo "    $(COMP_CC): "$<; \
	$(COMP_CC) $(CFLAGS_CC) $<
	@echo "#touch" > $(LIB_DIR)/lib.junk

.C.o:
	@ echo "    $(COMP_CC): "$<; \
	$(COMP_CC) $(CFLAGS_CC) $<
	@echo "#touch" > $(LIB_DIR)/lib.junk

.f.o:
	@ echo "    $(COMP_F): "$<; \
	$(COMP_F) $(CFLAGS_F) $<
	@echo "#touch" > $(LIB_DIR)/lib.junk

.F.o:
	@ echo "    $(COMP_F): "$<; \
	$(COMP_F) $(CFLAGS_F) $<
	@echo "#touch" > $(LIB_DIR)/lib.junk

.h.h_link:
	@ echo "    $(LN): "$<; \
	$(LN) $(CURR_DIR)/$< $(INC_DIR)/$<

.c.d:
	@ echo "    $(MAKE_DEPEND): "$<; \
	echo "" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.cpp.d:
	@ echo "    $(MAKE_DEPEND): "$<; \
	echo "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.C.d:
	@ echo "    $(MAKE_DEPEND): "$<; \
	echo "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.f.d:
	@ echo "    $(MAKE_DEPEND): "$<; \
	echo "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.F.d:
	@ echo "    $(MAKE_DEPEND): "$<; \
	echo "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

# dependencies
DEPEND = /dev/null
include $(DEPEND)
