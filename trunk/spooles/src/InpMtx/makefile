# $Id: makefile,v 1.2 2001-02-06 03:12:56 paklein Exp $
# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .C .f .F .o .d .h .h_link

# objects
OBJ = \
	InpMtx.adjForATA.o \
	InpMtx.basics.o \
	InpMtx.extract.o \
	InpMtx.fullAdj.o \
	InpMtx.gmmm.o \
	InpMtx.gmvm.o \
	InpMtx.init.o \
	InpMtx.input.o \
	InpMtx.instance.o \
	InpMtx.IO.o \
	InpMtx.map.o \
	InpMtx.mvm.o \
	InpMtx.mvmVector.o \
	InpMtx.permute.o \
	InpMtx.profile.o \
	InpMtx.support.o \
	InpMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	InpMtx.adjForATA.c \
	InpMtx.basics.c \
	InpMtx.extract.c \
	InpMtx.fullAdj.c \
	InpMtx.gmmm.c \
	InpMtx.gmvm.c \
	InpMtx.init.c \
	InpMtx.input.c \
	InpMtx.instance.c \
	InpMtx.IO.c \
	InpMtx.map.c \
	InpMtx.mvm.c \
	InpMtx.mvmVector.c \
	InpMtx.permute.c \
	InpMtx.profile.c \
	InpMtx.support.c \
	InpMtx.util.c

# headers
HDR = \
	InpMtx.h

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
	"DEPEND = $(DEPEND)" \
	"ECHO = $(ECHO)" \
	"CAT = $(CAT)"

# debugging
echo_all:
	@ $(ECHO) "MAKE = $(MAKE)"
	@ $(ECHO) "MAKE_OPTS = $(MAKE_OPTS)"
	@ $(ECHO) "INC_DIR = $(INC_DIR)"
	@ $(ECHO) "OBJ_DIR = $(OBJ_DIR)"
	@ $(ECHO) "LIB_DIR = $(LIB_DIR)"
	@ $(ECHO) "TR = $(TR)"
	@ $(ECHO) "COMP_C = $(COMP_C)"
	@ $(ECHO) "COMP_CC = $(COMP_CC)"
	@ $(ECHO) "COMP_F = $(COMP_F)"
	@ $(ECHO) "CFLAGS_C = $(CFLAGS_C)"
	@ $(ECHO) "CFLAGS_CC = $(CFLAGS_CC)"
	@ $(ECHO) "CFLAGS_F = $(CFLAGS_F)"
	@ $(ECHO) "LIBRARY = $(LIBRARY)"
	@ $(ECHO) "AR = $(AR)"
	@ $(ECHO) "ARFLAGS = $(ARFLAGS)"
	@ $(ECHO) "RM = $(RM)"
	@ $(ECHO) "RM_FILES = $(RM_FILES)"
	@ $(ECHO) "LN = $(LN)"
	@ $(ECHO) "MAKE_DEPEND = $(MAKE_DEPEND)"
	@ $(ECHO) "DEPEND_PATH = $(DEPEND_PATH)"
	@ $(ECHO) "DEPEND = $(DEPEND)"
	@ $(ECHO) "ECHO = $(ECHO)"
	@ $(ECHO) "CAT = $(CAT)"

clean:
	-@ $(RM) $(RM_FILES)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = clean"

depend_init: all.depend
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = depend_init"

all.depend: $(DEP)
	@ $(ECHO) "    $(CAT): $(DEP) > all.depend"; $(CAT) $(DEP) > all.depend

headers: $(HDR_LINK)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = headers"

objects: all.depend
	@ $(MAKE) $(MACROS) update_objects
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = objects"

update_objects: $(OBJ)

library:
	@ $(ECHO) "\\" >> $(LIBRARY)
	@ $(ECHO) "\t$(OBJ_LINK)\c" >> $(LIBRARY)
	@ $(MAKE) $(MACROS) subdir_driver "COMMAND = library"

$(SUB_DIR): dummy
	@ $(ECHO) " cd $@"; \
	cd $@; $(MAKE) $(COMMAND) $(MACROS) "CURR_DIR = $(CURR_DIR)/$@"

dummy:

# suffix rules

.c.o:
	@ $(ECHO) "    $(COMP_C): "$<; \
	$(COMP_C) $(CFLAGS_C) $<
	@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@
	@$(ECHO) "#touch" > $(LIB_DIR)/lib.junk

.cpp.o:
	@ $(ECHO) "    $(COMP_CC): "$<; \
	$(COMP_CC) $(CFLAGS_CC) $<
	@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@
	@$(ECHO) "#touch" > $(LIB_DIR)/lib.junk

.C.o:
	@ $(ECHO) "    $(COMP_CC): "$<; \
	$(COMP_CC) $(CFLAGS_CC) $<
	@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@
	@$(ECHO) "#touch" > $(LIB_DIR)/lib.junk

.f.o:
	@ $(ECHO) "    $(COMP_F): "$<; \
	$(COMP_F) $(CFLAGS_F) $<
	@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@
	@$(ECHO) "#touch" > $(LIB_DIR)/lib.junk

.F.o:
	@ $(ECHO) "    $(COMP_F): "$<; \
	$(COMP_F) $(CFLAGS_F) $<
	@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@
	@$(ECHO) "#touch" > $(LIB_DIR)/lib.junk

.h.h_link:
	@ $(ECHO) "    $(LN): "$<; \
	$(LN) $(CURR_DIR)/$< $(INC_DIR)/$<

.c.d:
	@ $(ECHO) "    $(MAKE_DEPEND): "$<; \
	$(ECHO) "" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.cpp.d:
	@ $(ECHO) "    $(MAKE_DEPEND): "$<; \
	$(ECHO) "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.C.d:
	@ $(ECHO) "    $(MAKE_DEPEND): "$<; \
	$(ECHO) "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.f.d:
	@ $(ECHO) "    $(MAKE_DEPEND): "$<; \
	$(ECHO) "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

.F.d:
	@ $(ECHO) "    $(MAKE_DEPEND): "$<; \
	$(ECHO) "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<

# dependencies
DEPEND = /dev/null
include $(DEPEND)
