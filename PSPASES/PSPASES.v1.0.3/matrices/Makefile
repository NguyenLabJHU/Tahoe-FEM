CC = cc

# for SGI, Cray
#FC = f90
# for IBM
FC = xlf

S1 = fcc2bin.c
T1 = fcc2bin
S2 = rsa2fcc.f
T2 = rsa2fcc

all : $(T1) $(T2)

$(T1) : $(S1)
	$(CC) -o $@ $(S1)

$(T2) : $(S2)
	$(FC) -o $@ $(S2)

clean :
	rm -f $(T1) $(T2)
