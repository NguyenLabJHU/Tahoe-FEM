# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	StringT.o \
	dMatrixT.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	StringT.cpp \
	dMatrixT.cpp

# headers
HDR = \
	StringT.h \
	dMatrixT.h \
	pArrayT.h \
	AutoFill2DT.h \
	nMatrixT.h \
	iosfwd.h \
	ios_fwd_decl.h \
	Constants.h \
	Environment.h \
	ExceptionCodes.h

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

StringT.o: /usr/include/math.h /usr/include/sgidefs.h
StringT.o: /usr/include/standards.h /usr/include/stdlib.h
StringT.o: /usr/include/string.h /usr/include/ctype.h
StringT.o: /usr/include/locale_attr.h /usr/include/CC/iostream.h
StringT.o: /usr/include/memory.h /usr/include/stddef.h
StringT.o: /usr/include/CC/fstream.h /usr/include/CC/iomanip.h
StringT.o: /usr/include/CC/generic.h Constants.h ExceptionCodes.h StringT.h
StringT.o: Environment.h /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
StringT.o: ios_fwd_decl.h
dMatrixT.o: /usr/include/CC/iostream.h /usr/include/memory.h
dMatrixT.o: /usr/include/stddef.h /usr/include/sgidefs.h
dMatrixT.o: /usr/include/CC/iomanip.h /usr/include/CC/generic.h Constants.h
dMatrixT.o: dMatrixT.h nMatrixT.h
dMatrixT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/nArrayT.h
dMatrixT.o: /usr/include/stdlib.h /usr/include/standards.h
dMatrixT.o: /usr/include/math.h
dMatrixT.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h Environment.h
dMatrixT.o: /usr/include/string.h ExceptionCodes.h
