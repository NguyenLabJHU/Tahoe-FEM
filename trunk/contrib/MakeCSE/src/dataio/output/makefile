# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	ExodusOutputT.o \
	FE_ASCIIT.o \
	OutputBaseT.o \
	OutputSetT.o \
	EnSightOutputT.o \
	AbaqusOutputT.o \
	TecPlotOutputT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	ExodusOutputT.cpp \
	FE_ASCIIT.cpp \
	OutputBaseT.cpp \
	OutputSetT.cpp \
	EnSightOutputT.cpp \
	AbaqusOutputT.cpp \
	TecPlotOutputT.cpp

# headers
HDR = \
	ExodusOutputT.h \
	FE_ASCIIT.h \
	OutputBaseT.h \
	OutputSetT.h \
	EnSightOutputT.h \
	AbaqusOutputT.h \
	TecPlotOutputT.h

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

ExodusOutputT.o: ExodusOutputT.h OutputBaseT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
ExodusOutputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
ExodusOutputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
ExodusOutputT.o: /usr/include/string.h /usr/include/standards.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
ExodusOutputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
ExodusOutputT.o: /usr/include/stdlib.h /usr/include/math.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusT.h
ExodusOutputT.o: OutputSetT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
ExodusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
ExodusOutputT.o: /usr/include/time.h
FE_ASCIIT.o: FE_ASCIIT.h OutputBaseT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
FE_ASCIIT.o: /usr/include/CC/iostream.h /usr/include/memory.h
FE_ASCIIT.o: /usr/include/stddef.h /usr/include/sgidefs.h
FE_ASCIIT.o: /usr/include/string.h /usr/include/standards.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
FE_ASCIIT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
FE_ASCIIT.o: /usr/include/stdlib.h /usr/include/math.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
FE_ASCIIT.o: /usr/include/CC/fstream.h OutputSetT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ModelFileT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
FE_ASCIIT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
OutputBaseT.o: OutputBaseT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
OutputBaseT.o: /usr/include/CC/iostream.h /usr/include/memory.h
OutputBaseT.o: /usr/include/stddef.h /usr/include/sgidefs.h
OutputBaseT.o: /usr/include/string.h /usr/include/standards.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
OutputBaseT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
OutputBaseT.o: /usr/include/stdlib.h /usr/include/math.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
OutputBaseT.o: OutputSetT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ModelFileT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/TriT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryBaseT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/QuadT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/HexahedronT.h
OutputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/TetrahedronT.h
OutputSetT.o: OutputSetT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
OutputSetT.o: /usr/include/CC/iostream.h /usr/include/memory.h
OutputSetT.o: /usr/include/stddef.h /usr/include/sgidefs.h
OutputSetT.o: /usr/include/string.h /usr/include/standards.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
OutputSetT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
OutputSetT.o: /usr/include/stdlib.h /usr/include/math.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
OutputSetT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
EnSightOutputT.o: EnSightOutputT.h OutputBaseT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
EnSightOutputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
EnSightOutputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
EnSightOutputT.o: /usr/include/string.h /usr/include/standards.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
EnSightOutputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
EnSightOutputT.o: /usr/include/stdlib.h /usr/include/math.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/EnSightT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
EnSightOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
EnSightOutputT.o: OutputSetT.h /usr/include/CC/fstream.h
AbaqusOutputT.o: AbaqusOutputT.h OutputBaseT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
AbaqusOutputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
AbaqusOutputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
AbaqusOutputT.o: /usr/include/string.h /usr/include/standards.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
AbaqusOutputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
AbaqusOutputT.o: /usr/include/stdlib.h /usr/include/math.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AbaqusT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h OutputSetT.h
AbaqusOutputT.o: /usr/include/CC/fstream.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
AbaqusOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
TecPlotOutputT.o: TecPlotOutputT.h OutputBaseT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
TecPlotOutputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
TecPlotOutputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
TecPlotOutputT.o: /usr/include/string.h /usr/include/standards.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
TecPlotOutputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
TecPlotOutputT.o: /usr/include/stdlib.h /usr/include/math.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/TecPlotT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
TecPlotOutputT.o: OutputSetT.h /usr/include/CC/fstream.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
TecPlotOutputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
