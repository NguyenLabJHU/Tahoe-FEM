# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	ExodusT.o \
	ModelFileT.o \
	EnSightT.o \
	AbaqusT.o \
	TecPlotT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	ExodusT.cpp \
	ModelFileT.cpp \
	EnSightT.cpp \
	AbaqusT.cpp \
	TecPlotT.cpp

# headers
HDR = \
	ExodusT.h \
	ModelFileT.h \
	EnSightT.h \
	AbaqusT.h \
	TecPlotT.h

# subdirectories
#SUBDIR = 

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

depend:
	-@ makedepend $(DEPENDPATH) $(SRC)

headers:
	-@ cp $(HDR) $(INCDIR)

update_headers: $(HEADERTEST)

touch:
	-@ touch -c $(SRC)
	-@ touch -c $(HDR)

# object file sets
objects: $(OBJ)

objects_1: $(OBJ_1)

objects_2: $(OBJ_2)

objects_3: $(OBJ_3)


$(HEADERTEST): $(HDR)
	@ echo "    cp: "$?; cp $? $(INCDIR)
	@ touch $(HEADERTEST)

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

ExodusT.o: ExodusT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
ExodusT.o: /usr/include/CC/iostream.h /usr/include/memory.h
ExodusT.o: /usr/include/stddef.h /usr/include/sgidefs.h
ExodusT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
ExodusT.o: /usr/include/stdlib.h /usr/include/standards.h /usr/include/math.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
ExodusT.o: /usr/include/string.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
ExodusT.o: /usr/include/time.h /supp/Access.new/inc/exodusII.h
ExodusT.o: /supp/Access.new/inc/netcdf.h /usr/include/errno.h
ExodusT.o: /usr/include/sys/errno.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/pArrayT.h
ExodusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ifstream_x.h
ExodusT.o: /usr/include/CC/fstream.h
ModelFileT.o: ModelFileT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
ModelFileT.o: /usr/include/CC/iostream.h /usr/include/memory.h
ModelFileT.o: /usr/include/stddef.h /usr/include/sgidefs.h
ModelFileT.o: /usr/include/string.h /usr/include/standards.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
ModelFileT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
ModelFileT.o: /usr/include/stdlib.h /usr/include/math.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
ModelFileT.o: /usr/include/ctype.h /usr/include/locale_attr.h
ModelFileT.o: /usr/include/float.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ifstream_x.h
ModelFileT.o: /usr/include/CC/fstream.h ExodusT.h
ModelFileT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
EnSightT.o: EnSightT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
EnSightT.o: /usr/include/CC/iostream.h /usr/include/memory.h
EnSightT.o: /usr/include/stddef.h /usr/include/sgidefs.h
EnSightT.o: /usr/include/string.h /usr/include/standards.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
EnSightT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
EnSightT.o: /usr/include/stdlib.h /usr/include/math.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ifstream_x.h
EnSightT.o: /usr/include/CC/fstream.h /usr/include/ctype.h
EnSightT.o: /usr/include/locale_attr.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
EnSightT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
AbaqusT.o: AbaqusT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
AbaqusT.o: /usr/include/CC/iostream.h /usr/include/memory.h
AbaqusT.o: /usr/include/stddef.h /usr/include/sgidefs.h /usr/include/string.h
AbaqusT.o: /usr/include/standards.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
AbaqusT.o: /usr/include/stdlib.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
AbaqusT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
AbaqusT.o: /usr/include/math.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
AbaqusT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
AbaqusT.o: /usr/include/time.h
TecPlotT.o: TecPlotT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
TecPlotT.o: /usr/include/CC/iostream.h /usr/include/memory.h
TecPlotT.o: /usr/include/stddef.h /usr/include/sgidefs.h
TecPlotT.o: /usr/include/string.h /usr/include/standards.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
TecPlotT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
TecPlotT.o: /usr/include/stdlib.h /usr/include/math.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
TecPlotT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
