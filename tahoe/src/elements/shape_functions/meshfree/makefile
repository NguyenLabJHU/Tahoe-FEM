# $Id: makefile,v 1.3 2001-02-06 03:13:13 paklein Exp $
# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .C .f .F .o .d .h .h_link

# objects
OBJ = \
	MeshFreeShapeFunctionT.o \
	MeshFreeSupportT.o \
	MeshFreeSupport2DT.o \
	MeshFreeSupport3DT.o \
	MeshFreeSurfaceShapeT.o \
	MeshFreeSurfaceSupportT.o \
	MeshFreeT.o \
	MLSSolverT.o \
	OrthoMLSSolverT.o \
	OrthoMLS2DT.o \
	OrthoMLS3DT.o \
	SamplingSurfaceT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MeshFreeShapeFunctionT.cpp \
	MeshFreeSupportT.cpp \
	MeshFreeSupport2DT.cpp \
	MeshFreeSupport3DT.cpp \
	MeshFreeSurfaceShapeT.cpp \
	MeshFreeSurfaceSupportT.cpp \
	MeshFreeT.cpp \
	MLSSolverT.cpp \
	OrthoMLSSolverT.cpp \
	OrthoMLS2DT.cpp \
	OrthoMLS3DT.cpp \
	SamplingSurfaceT.cpp

# headers
HDR = \
	MeshFreeShapeFunctionT.h \
	MeshFreeSupportT.h \
	MeshFreeSupport2DT.h \
	MeshFreeSupport3DT.h \
	MeshFreeSurfaceShapeT.h \
	MeshFreeSurfaceSupportT.h \
	MeshFreeT.h \
	MLSSolverT.h \
	OrthoMLSSolverT.h \
	OrthoMLS2DT.h \
	OrthoMLS3DT.h \
	SamplingSurfaceT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	basis \
	D2
#SUB_DIR = no_subdirs

# process sub-directories as targets:
subdir_driver: $(SUB_DIR)
#subdir_driver:

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
