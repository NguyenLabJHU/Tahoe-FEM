# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	GeometryBaseT.o \
	GeometryT.o \
	HexahedronT.o \
	LineT.o \
	PentahedronT.o \
	QuadT.o \
	TetrahedronT.o \
	TriT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	GeometryBaseT.cpp \
	GeometryT.cpp \
	HexahedronT.cpp \
	LineT.cpp \
	PentahedronT.cpp \
	QuadT.cpp \
	TetrahedronT.cpp \
	TriT.cpp

# headers
HDR = \
	GeometryBaseT.h \
	GeometryT.h \
	HexahedronT.h \
	LineT.h \
	PentahedronT.h \
	QuadT.h \
	TetrahedronT.h \
	TriT.h

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

GeometryBaseT.o: GeometryBaseT.h GeometryT.h
GeometryBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
GeometryBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
GeometryT.o: GeometryT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
GeometryT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
GeometryT.o: /usr/include/CC/iostream.h /usr/include/memory.h
GeometryT.o: /usr/include/stddef.h /usr/include/sgidefs.h
GeometryT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
HexahedronT.o: /usr/include/math.h /usr/include/sgidefs.h
HexahedronT.o: /usr/include/standards.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
HexahedronT.o: HexahedronT.h GeometryBaseT.h GeometryT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
HexahedronT.o: /usr/include/CC/iostream.h /usr/include/memory.h
HexahedronT.o: /usr/include/stddef.h /usr/include/CC/iomanip.h
HexahedronT.o: /usr/include/CC/generic.h /usr/include/stdlib.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
HexahedronT.o: /usr/include/string.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
HexahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
LineT.o: /usr/include/math.h /usr/include/sgidefs.h /usr/include/standards.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h LineT.h
LineT.o: GeometryBaseT.h GeometryT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
LineT.o: /usr/include/CC/iostream.h /usr/include/memory.h
LineT.o: /usr/include/stddef.h /usr/include/CC/iomanip.h
LineT.o: /usr/include/CC/generic.h /usr/include/stdlib.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
LineT.o: /usr/include/string.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
LineT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
PentahedronT.o: PentahedronT.h GeometryBaseT.h GeometryT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h QuadT.h
PentahedronT.o: TriT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
PentahedronT.o: /usr/include/CC/iostream.h /usr/include/memory.h
PentahedronT.o: /usr/include/stddef.h /usr/include/sgidefs.h
PentahedronT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
PentahedronT.o: /usr/include/stdlib.h /usr/include/standards.h
PentahedronT.o: /usr/include/math.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
PentahedronT.o: /usr/include/string.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
PentahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
QuadT.o: /usr/include/math.h /usr/include/sgidefs.h /usr/include/standards.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h QuadT.h
QuadT.o: GeometryBaseT.h GeometryT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
QuadT.o: /usr/include/CC/iostream.h /usr/include/memory.h
QuadT.o: /usr/include/stddef.h /usr/include/CC/iomanip.h
QuadT.o: /usr/include/CC/generic.h /usr/include/stdlib.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
QuadT.o: /usr/include/string.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
QuadT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
TetrahedronT.o: TetrahedronT.h GeometryBaseT.h GeometryT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h QuadT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
TetrahedronT.o: /usr/include/CC/iostream.h /usr/include/memory.h
TetrahedronT.o: /usr/include/stddef.h /usr/include/sgidefs.h
TetrahedronT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
TetrahedronT.o: /usr/include/stdlib.h /usr/include/standards.h
TetrahedronT.o: /usr/include/math.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
TetrahedronT.o: /usr/include/string.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
TetrahedronT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h TriT.h
TriT.o: GeometryBaseT.h GeometryT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h QuadT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
TriT.o: /usr/include/CC/iostream.h /usr/include/memory.h
TriT.o: /usr/include/stddef.h /usr/include/sgidefs.h
TriT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
TriT.o: /usr/include/stdlib.h /usr/include/standards.h /usr/include/math.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h /usr/include/string.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dMatrixT.h
TriT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nMatrixT.h
