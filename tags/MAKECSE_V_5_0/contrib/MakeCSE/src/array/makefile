# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	dArrayT.o \
	dArray2DT.o \
	iArrayT.o \
	iArray2DT.o \
	iAutoArrayT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	dArrayT.cpp \
	dArray2DT.cpp \
	iArrayT.cpp \
	iArray2DT.cpp \
	iAutoArrayT.cpp

# headers
HDR = \
	ArrayT.h \
	Array2DT.h \
	AutoArrayT.h \
	dArrayT.h \
	dArray2DT.h \
	iArrayT.h \
	iArray2DT.h \
	iAutoArrayT.h \
	nArrayT.h \
	nArray2DT.h \
	RaggedArray2DT.h \
	VariRaggedArray2DT.h

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

dArrayT.o: /usr/include/CC/iostream.h /usr/include/memory.h
dArrayT.o: /usr/include/stddef.h /usr/include/sgidefs.h
dArrayT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
dArrayT.o: /usr/include/math.h /usr/include/standards.h
dArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h dArrayT.h
dArrayT.o: nArrayT.h /usr/include/stdlib.h ArrayT.h
dArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
dArrayT.o: /usr/include/string.h
dArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
dArray2DT.o: /usr/include/CC/iostream.h /usr/include/memory.h
dArray2DT.o: /usr/include/stddef.h /usr/include/sgidefs.h
dArray2DT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
dArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h dArray2DT.h
dArray2DT.o: nArray2DT.h nArrayT.h /usr/include/stdlib.h
dArray2DT.o: /usr/include/standards.h /usr/include/math.h ArrayT.h
dArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
dArray2DT.o: /usr/include/string.h
dArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
dArray2DT.o: iArrayT.h dArrayT.h
iArrayT.o: /usr/include/CC/iostream.h /usr/include/memory.h
iArrayT.o: /usr/include/stddef.h /usr/include/sgidefs.h
iArrayT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
iArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h iArrayT.h
iArrayT.o: nArrayT.h /usr/include/stdlib.h /usr/include/standards.h
iArrayT.o: /usr/include/math.h ArrayT.h
iArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
iArrayT.o: /usr/include/string.h
iArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
iArray2DT.o: /usr/include/CC/iostream.h /usr/include/memory.h
iArray2DT.o: /usr/include/stddef.h /usr/include/sgidefs.h
iArray2DT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
iArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h iArray2DT.h
iArray2DT.o: nArray2DT.h nArrayT.h /usr/include/stdlib.h
iArray2DT.o: /usr/include/standards.h /usr/include/math.h ArrayT.h
iArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
iArray2DT.o: /usr/include/string.h
iArray2DT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
iAutoArrayT.o: iAutoArrayT.h AutoArrayT.h ArrayT.h
iAutoArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
iAutoArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
iAutoArrayT.o: /usr/include/CC/iostream.h /usr/include/memory.h
iAutoArrayT.o: /usr/include/stddef.h /usr/include/sgidefs.h
iAutoArrayT.o: /usr/include/string.h /usr/include/standards.h
iAutoArrayT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
iAutoArrayT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h
