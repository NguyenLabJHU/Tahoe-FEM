include ../Makefile.in


CFLAGS = $(PCFLAGS) -I. $(INCDIR)


OBJS = comm.o util.o memory.o debug.o setup.o \
       kmetis.o rmetis.o diffuse.o gmetis.o ometis.o scremap.o grsetup.o timer.o \
       coarsen.o edge_refine.o initpart.o drivers.o fold.o \
       initdiff.o wave.o move.o node_refine.o initmsection.o order.o pqueue.o \
       diffutil.o xyzpart.o pspases.o remap.o frename.o \
       iintsort.o iidxsort.o ikeysort.o ikeyvalsort.o

.c.o:
	$(PCC) $(CFLAGS) -c $*.c


libparmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o libparmetis.a
