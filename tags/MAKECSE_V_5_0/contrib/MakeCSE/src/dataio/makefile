# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	IOBaseT.o \
	IOManager.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	IOBaseT.cpp \
	IOManager.cpp 

# headers
HDR = \
	IOBaseT.h \
	IOManager.h \
	MapSetT.h

# subdirectories
SUBDIR = \
	database \
	input \
	output

MACROS = \
	"MAKEOPTS    = $(MAKEOPTS)" \
	"INCDIR      = $(INCDIR)" \
	"OBJDIR      = $(OBJDIR)" \
	"TR          = $(TR)" \
	"COMP_CC     = $(COMP_CC)" \
	"COMP_C      = $(COMP_C)" \
	"CFLAGS      = $(CFLAGS)" \
	"MORE_CFLAGS = $(MORE_CFLAGS)" \
	"DEPENDPATH  = $(DEPENDPATH)"

HEADERTEST = makefile

# drive subdirectories as targets:
subdir_driver: $(SUBDIR)

clean:
	-@ rm -r *.o *.i *.*~ *.bak ii* *.o1 *.o2 *.o3
	@ make $(MAKEOPTS) subdir_driver "COMMAND = clean"

depend:
	-@ makedepend $(DEPENDPATH) $(SRC)
	@ make $(MAKEOPTS) subdir_driver "COMMAND = depend"

headers:
	-@ cp $(HDR) $(INCDIR)
	@ make $(MAKEOPT) subdir_driver "COMMAND = headers"

update_headers: $(HEADERTEST)
	@ make $(MAKEOPT) subdir_driver "COMMAND = update_headers"

touch:
	-@ touch -c $(SRC)
	-@ touch -c $(HDR)
	@ make $(MAKEOPT) subdir_driver "COMMAND = touch"

# object file sets
objects: $(OBJ)
	@ make $(MAKEOPTS) subdir_driver "COMMAND = objects"
objects_1: $(OBJ_1)
	@ make $(MAKEOPTS) subdir_driver "COMMAND = objects_1"
objects_2: $(OBJ_2)
	@ make $(MAKEOPTS) subdir_driver "COMMAND = objects_2"
objects_3: $(OBJ_3)
	@ make $(MAKEOPTS) subdir_driver "COMMAND = objects_3"

$(HEADERTEST): $(HDR)
	@ echo "    cp: "$?; cp $? $(INCDIR)
	@ touch $(HEADERTEST)

$(SUBDIR): dummy
	@ echo " cd $@"; \
	cd $@; make $(COMMAND) $(MACROS)

dummy:

# suffix rules
DO_COMP_CC = \
	@ echo "    $(COMP_CC): "$<; \
	$(COMP_CC) $(CFLAGS) $(MORE_CFLAGS) $<; \
	cp $@ $(OBJDIR); \
	touch $(OBJDIR)/obj.junk

.cpp.o:
	$(DO_COMP_CC)
.cpp.o1:
	$(DO_COMP_CC)
.cpp.o2:
	$(DO_COMP_CC)
.cpp.o3:
	$(DO_COMP_CC)

DO_COMP_C = \
	@ echo "    $(COMP_C): "$<; \
	$(COMP_C) $(CFLAGS) $(MORE_CFLAGS) $<; \
	cp $@ $(OBJDIR); \
	touch $(OBJDIR)/obj.junk

.c.o:
	$(DO_COMP_C)
.c.o1:
	$(DO_COMP_C)
.c.o2:
	$(DO_COMP_C)
.c.o3:
	$(DO_COMP_C)

# DO NOT DELETE

IOBaseT.o: IOBaseT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
IOBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
IOBaseT.o: /usr/include/CC/iostream.h /usr/include/memory.h
IOBaseT.o: /usr/include/stddef.h /usr/include/sgidefs.h
IOBaseT.o: /usr/include/CC/fstream.h /usr/include/CC/iomanip.h
IOBaseT.o: /usr/include/CC/generic.h
IOBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
IOManager.o: IOManager.h /usr/include/CC/iostream.h /usr/include/memory.h
IOManager.o: /usr/include/stddef.h /usr/include/sgidefs.h
IOManager.o: /usr/include/CC/fstream.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
IOManager.o: /usr/include/string.h /usr/include/standards.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
IOManager.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
IOManager.o: /usr/include/stdlib.h /usr/include/math.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h IOBaseT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ifstream_x.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/OutputSetT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/InputBaseT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusInputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/TextInput.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ModelFileT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/EnSightInputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/EnSightT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AbaqusInputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AbaqusT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/FE_ASCIIT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/OutputBaseT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusOutputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/EnSightOutputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AbaqusOutputT.h
IOManager.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/TecPlotOutputT.h
