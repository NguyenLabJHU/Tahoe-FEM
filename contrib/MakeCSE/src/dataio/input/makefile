# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	ExodusInputT.o \
	InputBaseT.o \
	TextInput.o \
	EnSightInputT.o \
	AbaqusInputT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	ExodusInputT.cpp \
	InputBaseT.cpp \
	TextInput.cpp \
	EnSightInputT.cpp \
	AbaqusInputT.cpp

# headers
HDR = \

	ExodusInputT.h \
	InputBaseT.h \
	TextInput.h \
	EnSightInputT.h \
	AbaqusInputT.h

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

ExodusInputT.o: InputBaseT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
ExodusInputT.o: ExodusInputT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExodusT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
ExodusInputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
ExodusInputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
ExodusInputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
ExodusInputT.o: /usr/include/stdlib.h /usr/include/standards.h
ExodusInputT.o: /usr/include/math.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
ExodusInputT.o: /usr/include/string.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
ExodusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
InputBaseT.o: InputBaseT.h /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
InputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
InputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
InputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
InputBaseT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
TextInput.o: TextInput.h InputBaseT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ModelFileT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
TextInput.o: /usr/include/CC/iostream.h /usr/include/memory.h
TextInput.o: /usr/include/stddef.h /usr/include/sgidefs.h
TextInput.o: /usr/include/string.h /usr/include/standards.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
TextInput.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
TextInput.o: /usr/include/stdlib.h /usr/include/math.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
TextInput.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
EnSightInputT.o: EnSightInputT.h InputBaseT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/EnSightT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
EnSightInputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
EnSightInputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
EnSightInputT.o: /usr/include/string.h /usr/include/standards.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
EnSightInputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
EnSightInputT.o: /usr/include/stdlib.h /usr/include/math.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ifstream_x.h
EnSightInputT.o: /usr/include/CC/fstream.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
EnSightInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
AbaqusInputT.o: AbaqusInputT.h InputBaseT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/IOBaseT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/GeometryT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AbaqusT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iosfwd.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
AbaqusInputT.o: /usr/include/CC/iostream.h /usr/include/memory.h
AbaqusInputT.o: /usr/include/stddef.h /usr/include/sgidefs.h
AbaqusInputT.o: /usr/include/string.h /usr/include/standards.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArray2DT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArray2DT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
AbaqusInputT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
AbaqusInputT.o: /usr/include/stdlib.h /usr/include/math.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iArrayT.h
AbaqusInputT.o: /usr/include/CC/fstream.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArray2DT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/dArrayT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/iAutoArrayT.h
AbaqusInputT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/AutoArrayT.h
