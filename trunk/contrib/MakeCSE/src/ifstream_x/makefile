# object directory makefile - compilation instructions

.SUFFIXES: .c .cpp .o .o1 .o2 .o3

# objects
OBJ = \
	ifstream_x.o

OBJ_1 = $(OBJ:.o=.o1)
OBJ_2 = $(OBJ:.o=.o2)
OBJ_3 = $(OBJ:.o=.o3)

# sources
SRC = \
	ifstream_x.cpp

# headers
HDR = \
	ifstream_x.h

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

ifstream_x.o: ifstream_x.h /usr/include/CC/fstream.h
ifstream_x.o: /usr/include/CC/iostream.h /usr/include/memory.h
ifstream_x.o: /usr/include/stddef.h /usr/include/sgidefs.h
ifstream_x.o: /usr/include/string.h /usr/include/standards.h
ifstream_x.o: /usr/include/ctype.h /usr/include/locale_attr.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Environment.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ExceptionCodes.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/StringT.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ArrayT.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/Constants.h
ifstream_x.o: /u/jwfoulk/sawimme_fracture/tetexo/inc/ios_fwd_decl.h
