.SUFFIXES: .c .f

F_FILES = alltoallu_hc.f \
alltoonev.f \
b_ax.f \
dpotrf.f \
emovea.f \
eparfact1i.f \
gen_lc.f \
genleaftree.f \
gentree.f \
getmyhvb.f \
getmysnodes.f \
moveai.f \
moveav.f \
order.f \
parelimh1i.f \
parelimv1i.f \
parfact1i.f \
parsymb.f \
pbsolve1.f \
pbsolvem.f \
pfsolve1.f \
pfsolvem.f \
pmovea.f \
pparfact1i.f \
premovea.f \
preordbc.f \
preordbe.f \
preordxc.f \
reordb.f \
reordx.f \
serialfactor.f \
subpspaces.f \
trisolve.f \
ygentree.f

C_FILES = $(F_FILES:.f=.c)

F2C_OPTS = -I/usr/local/lam/include -c -w -a

all: $(C_FILES)

clean:
	rm $(C_FILES)

.f.c:
	f2c $(F2C_OPTS) $<
