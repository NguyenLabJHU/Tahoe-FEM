# --------------------- Compilers and Linkers -------------------------
## for MPICH (Useful for network of workstations, which do not have
##            special compiler support for MPI libraries)
CC = mpicc
FC = mpif90

## for IBM
#CC = mpcc
#FC = mpxlf 

## for Cray and Origin
#CC = cc 
#FC = f90 

## for 64-bit code on SGI -- tested okay on Origin 2000.
#CC = cc -64
#FC = cc -64

FLINKER = $(FC) 
CLINKER = $(CC) 

# --------------------- BLAS Selection        -------------------------

## for Intel's Math Kernel Library (MKL)
BLAS = -L/apps/blas -lblas
#BLAS = /apps/blas/mkl_libblas.a
#BLAS = -L/apps/intel/mkl/lib/32 -lmkl_lapack -lmkl_ia32 lib/libaux32_intel.a
#BLAS = -L/apps/intel/mkl/lib/32 -lmkl -lmkl_lapack32 -lmkl_def

## for native SGI BLAS
#BLAS = -lblas

## for native IBM BLAS (use esslp2 for POWER2 architecture. if thread-safe 
##                      compilers are used then esslsmp or essl_r can be used.)
# BLAS = -lessl

## for Cray, SciLib is automatically linked. if not, uncomment the following.
# BLAS = -lsci

## for SUN native BLAS (Using licensed Sun Performance Library)
# BLAS = -xlic_lib=sunperf 

## To supply your own BLAS, change BLASDIR value below by the directory where 
## you have built a BLAS library. Following assumes that you have built 
## libmyblas.a in directory ../../BLAS (relative to current, TEST, directory).
# BLASDIR = ../../BLAS
# BLAS = -L $(BLASDIR) -lmyblas

# --------------------- Driver Related Selection ----------------------

## for IBM
# FLFL = 
# CLFL = -lxlf90 

## for Cray
# PSTOBJS = craydep.o
# CFL = -DINT8

## for SGI (MIPSPro 7.2 F90 compilers)
# CLFL = -lftn -lfortran

## for SUN 
# CLFL = -lf90 -lF77 -lm

# --------------------- External libraries    -------------------------
## for SGI Origin 
# EXTLIB = -lmpi

# ---------------------------------------------------------------------

INCLUDE_DIR = -I.
PSPASESDIR = ../PSPASES
PARMETISDIR = ../parMetis
METISDIR = ../metis
LIB_DIR = -L$(METISDIR) -L$(PARMETISDIR) -L$(PSPASESDIR)
LIBS = -lpspases -lparmetis -lmetis -lm $(BLAS) $(EXTLIB)
LIB_DEP = $(PSPASESDIR)/libpspases.a $(PARMETISDIR)/libparmetis.a $(METISDIR)/libmetis.a

S1 = pspases_testfc.o $(PSTOBJS)
TARGET1 = fcpspases_test

S2 = pspases_testfr.o $(PSTOBJS)
TARGET2 = frpspases_test

S3 = pspases_testc.o 
TARGET3 = cpspases_test

all : $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1) : $(S1) $(LIB_DEP)
	$(FLINKER) -o $@ $(S1) $(LIB_DIR) $(LIBS) $(FLFL)

$(TARGET2) : $(S2) $(LIB_DEP)
	$(FLINKER) -o $@ $(S2) $(LIB_DIR) $(LIBS) $(FLFL)

$(TARGET3) : $(S3) $(LIB_DEP)
	$(CLINKER) -o $@ $(S3) $(LIB_DIR) $(LIBS) $(CLFL)

.f.o:
	$(FC) -c $(FFL) $(INCLUDE_DIR) $*.f
.c.o: 
	$(CC) -c $(CFL) $(INCLUDE_DIR) $*.c
clean:
	rm -f *.o $(TARGET1) $(TARGET2) $(TARGET3)

