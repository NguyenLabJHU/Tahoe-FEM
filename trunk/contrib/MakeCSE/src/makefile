# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	GlobalEdgeFinderT.o \
	ElementBaseT.o \
	CSEBaseT.o \
	FEManager.o \
	MakeCSE.o \
	MakeCSEIOManager.o \
	NodeManagerPrimitive.o \
	Quad2Tri.o \
	main.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	GlobalEdgeFinderT.cpp \
	ElementBaseT.cpp \
	CSEBaseT.cpp \
	FEManager.cpp \
	MakeCSE.cpp \
	MakeCSEIOManager.cpp \
	NodeManagerPrimitive.cpp \
	Quad2Tri.cpp \
	main.cpp

# headers
HDR = \
	GlobalEdgeFinderT.h \
	ElementBaseT.h \
	CSEBaseT.h \
	FEManager.h \
	MakeCSE.h \
	MakeCSEIOManager.h \
	NodeManagerPrimitive.h \
	Quad2Tri.h \
	ClockT.h

# subdirectories
SUBDIR = \
	array \
	misc \
	ifstream_x \
	geometry \
	dataio

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

GlobalEdgeFinderT.o: GlobalEdgeFinderT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
GlobalEdgeFinderT.o: /usr/include/CC/iostream.h /usr/include/memory.h
GlobalEdgeFinderT.o: /usr/include/stddef.h /usr/include/sgidefs.h
GlobalEdgeFinderT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
GlobalEdgeFinderT.o: /usr/include/stdlib.h /usr/include/standards.h
GlobalEdgeFinderT.o: /usr/include/math.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
GlobalEdgeFinderT.o: /usr/include/string.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
GlobalEdgeFinderT.o: FEManager.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
GlobalEdgeFinderT.o: MakeCSE.h ClockT.h /usr/include/time.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
GlobalEdgeFinderT.o: ElementBaseT.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
GlobalEdgeFinderT.o: NodeManagerPrimitive.h
GlobalEdgeFinderT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
ElementBaseT.o: ElementBaseT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
ElementBaseT.o: /usr/include/CC/iostream.h /usr/include/memory.h
ElementBaseT.o: /usr/include/stddef.h /usr/include/sgidefs.h
ElementBaseT.o: /usr/include/string.h /usr/include/standards.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
ElementBaseT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
ElementBaseT.o: /usr/include/stdlib.h /usr/include/math.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
ElementBaseT.o: FEManager.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h MakeCSE.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
ElementBaseT.o: ClockT.h /usr/include/time.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
ElementBaseT.o: GlobalEdgeFinderT.h MakeCSEIOManager.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
ElementBaseT.o: /usr/include/CC/fstream.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/TriT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryBaseT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/QuadT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/HexahedronT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/TetrahedronT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/PentahedronT.h
ElementBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/OutputSetT.h
CSEBaseT.o: CSEBaseT.h ElementBaseT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
CSEBaseT.o: /usr/include/CC/iostream.h /usr/include/memory.h
CSEBaseT.o: /usr/include/stddef.h /usr/include/sgidefs.h
CSEBaseT.o: /usr/include/string.h /usr/include/standards.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
CSEBaseT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
CSEBaseT.o: /usr/include/stdlib.h /usr/include/math.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/TriT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryBaseT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/QuadT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/HexahedronT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/TetrahedronT.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/PentahedronT.h
CSEBaseT.o: FEManager.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
CSEBaseT.o: MakeCSE.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
CSEBaseT.o: ClockT.h /usr/include/time.h
CSEBaseT.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
CSEBaseT.o: GlobalEdgeFinderT.h
FEManager.o: FEManager.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
FEManager.o: MakeCSE.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
FEManager.o: /usr/include/CC/iostream.h /usr/include/memory.h
FEManager.o: /usr/include/stddef.h /usr/include/sgidefs.h
FEManager.o: /usr/include/string.h /usr/include/standards.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
FEManager.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
FEManager.o: /usr/include/stdlib.h /usr/include/math.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h ClockT.h
FEManager.o: /usr/include/time.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
FEManager.o: GlobalEdgeFinderT.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
FEManager.o: MakeCSEIOManager.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
FEManager.o: /usr/include/CC/fstream.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
FEManager.o: NodeManagerPrimitive.h
FEManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
FEManager.o: ElementBaseT.h CSEBaseT.h Quad2Tri.h
MakeCSE.o: MakeCSE.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
MakeCSE.o: /usr/include/CC/iostream.h /usr/include/memory.h
MakeCSE.o: /usr/include/stddef.h /usr/include/sgidefs.h /usr/include/string.h
MakeCSE.o: /usr/include/standards.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
MakeCSE.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
MakeCSE.o: /usr/include/stdlib.h /usr/include/math.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h ClockT.h
MakeCSE.o: /usr/include/time.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
MakeCSE.o: FEManager.h GlobalEdgeFinderT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
MakeCSE.o: NodeManagerPrimitive.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
MakeCSE.o: ElementBaseT.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
MakeCSE.o: MakeCSEIOManager.h
MakeCSE.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
MakeCSE.o: /usr/include/CC/fstream.h
MakeCSEIOManager.o: MakeCSEIOManager.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
MakeCSEIOManager.o: /usr/include/CC/iostream.h /usr/include/memory.h
MakeCSEIOManager.o: /usr/include/stddef.h /usr/include/sgidefs.h
MakeCSEIOManager.o: /usr/include/CC/fstream.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
MakeCSEIOManager.o: /usr/include/string.h /usr/include/standards.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
MakeCSEIOManager.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
MakeCSEIOManager.o: /usr/include/stdlib.h /usr/include/math.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
MakeCSEIOManager.o: /usr/include/ctype.h /usr/include/locale_attr.h
MakeCSEIOManager.o: Quad2Tri.h ElementBaseT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
MakeCSEIOManager.o: FEManager.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
MakeCSEIOManager.o: MakeCSE.h ClockT.h /usr/include/time.h
MakeCSEIOManager.o: GlobalEdgeFinderT.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ifstream_x.h
MakeCSEIOManager.o: NodeManagerPrimitive.h
MakeCSEIOManager.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
NodeManagerPrimitive.o: NodeManagerPrimitive.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
NodeManagerPrimitive.o: /usr/include/CC/iostream.h /usr/include/memory.h
NodeManagerPrimitive.o: /usr/include/stddef.h /usr/include/sgidefs.h
NodeManagerPrimitive.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
NodeManagerPrimitive.o: /usr/include/stdlib.h /usr/include/standards.h
NodeManagerPrimitive.o: /usr/include/math.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
NodeManagerPrimitive.o: /usr/include/string.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
NodeManagerPrimitive.o: FEManager.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
NodeManagerPrimitive.o: MakeCSE.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
NodeManagerPrimitive.o: ClockT.h /usr/include/time.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
NodeManagerPrimitive.o: GlobalEdgeFinderT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
NodeManagerPrimitive.o: MakeCSEIOManager.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
NodeManagerPrimitive.o: /usr/include/CC/fstream.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
NodeManagerPrimitive.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArrayT.h
Quad2Tri.o: Quad2Tri.h ElementBaseT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
Quad2Tri.o: /usr/include/CC/iostream.h /usr/include/memory.h
Quad2Tri.o: /usr/include/stddef.h /usr/include/sgidefs.h
Quad2Tri.o: /usr/include/string.h /usr/include/standards.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
Quad2Tri.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
Quad2Tri.o: /usr/include/stdlib.h /usr/include/math.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArrayT.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/dArray2DT.h
Quad2Tri.o: NodeManagerPrimitive.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
Quad2Tri.o: FEManager.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
Quad2Tri.o: MakeCSE.h ClockT.h /usr/include/time.h
Quad2Tri.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
Quad2Tri.o: GlobalEdgeFinderT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iosfwd.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Environment.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/StringT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ArrayT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/Constants.h
main.o: /usr/include/CC/iostream.h /usr/include/memory.h
main.o: /usr/include/stddef.h /usr/include/sgidefs.h /usr/include/string.h
main.o: /usr/include/standards.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ExceptionCodes.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ios_fwd_decl.h FEManager.h
main.o: MakeCSE.h /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArray2DT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArray2DT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/nArrayT.h
main.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
main.o: /usr/include/stdlib.h /usr/include/math.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iArrayT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/iAutoArrayT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/AutoArrayT.h ClockT.h
main.o: /usr/include/time.h GlobalEdgeFinderT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOBaseT.h
main.o: MakeCSEIOManager.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/IOManager.h
main.o: /usr/include/CC/fstream.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/GeometryT.h
main.o: /u/jwfoulk/sawimme_fracture/MakeCSEv5/inc/ifstream_x.h
