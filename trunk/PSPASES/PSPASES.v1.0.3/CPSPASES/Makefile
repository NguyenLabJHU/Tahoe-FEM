include ../Makefile.in

FC=$(PFC)
CC=$(PCC)
FFLAGS=$(PFFLAGS)
CFLAGS=$(PCFLAGS)

INCLUDE_DIR = -I.

SOLVER_LIB = cpspases

SsolvF = trisolve.o order.o gentree.o genleaftree.o emovea.o pmovea.o parsymb.o alltoonev.o pbsolve1.o pbsolvem.o pfsolve1.o pfsolvem.o getmyhvb.o subpspaces.o parfact1i.o parelimh1i.o parelimv1i.o getmysnodes.o alltoallu_hc.o b_ax.o serialfactor.o gen_lc.o preordbe.o preordbc.o preordxc.o reordb.o reordx.o premovea.o moveai.o moveav.o eparfact1i.o pparfact1i.o ygentree.o psp_debug_f.o
SsolvC = util.o dpspacef.o dpspacet.o pspacec.o parmetisf.o auxil.o dpspacen.o pspaceo.o pspacey.o psp_debug.o
Ssolv = $(SsolvF) $(SsolvC) $(PSOBJS)

slib = lib$(SOLVER_LIB).a

all:	$(slib)

$(slib) : $(Ssolv)
	$(AR) $@ $(Ssolv) 

dpspacef.o:	pspaces.h
dpspacet.o:	pspaces.h
pspacec.o:	pspaces.h
auxil.o:	pspaces.h

dpotrf.o:
	$(FC) -c $(BLASFLAG) dpotrf.f

#$(FC) -c $(BLASFLAG) dpotrf.f

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIR) $*.f
.c.o: 
	$(CC) -c $(CFLAGS) $(INCLUDE_DIR) $*.c
clean:
	rm -f $(slib) $(Ssolv) $(BLASOBJ)
cleanc:
	rm -f $(SsolvC)
cleanf:
	rm -f $(SsolvF)
