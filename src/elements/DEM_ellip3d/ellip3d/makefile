# C++ compiler
CC = g++

# optimization (NDEBUG is only related to cassert)
OPTIMIZE = -O3 -DNDEBUG

# debugging options (macro DEBUG is defined in root6.cpp)
DEBUG = -Wall -g -DDEBUG

# OPENMP (macro OPENMP is defined in assembly.cpp )
OPENMP = -DOPENMP -fopenmp

# $(DEBUG) or $(OPTIMIZE), $(OPENMP) or not
OPTION = $(OPTIMIZE)
#OPTION = $(OPTIMIZE) $(OPENMP)

ellip3d: main.o assembly.o cylinder.o rectangle.o particle.o root6.o ran.o balanc.o zrhqr.o hqr.o vec.o nrutil.o parameter.o
	$(CC) $(OPTION) -o ellip3d main.o assembly.o cylinder.o rectangle.o particle.o root6.o ran.o balanc.o zrhqr.o hqr.o vec.o nrutil.o parameter.o

main.o: parameter.h gradation.h rectangle.h assembly.h main.cpp
	$(CC) -c $(OPTION) -o main.o main.cpp

assembly.o: assembly.h parameter.h vec.h gradation.h particle.h contact.h boundary.h rectangle.h cylinder.h assembly.cpp
	$(CC) -c $(OPTION) assembly.cpp

particle.o: particle.h parameter.h ran.h root6.h vec.h gradation.h contact.h boundary.h rectangle.h cylinder.h boundarytgt.h particle.cpp
	$(CC) -c $(OPTION) particle.cpp

cylinder.o: cylinder.h parameter.h ran.h vec.h shape.h cylinder.cpp
	$(CC) -c $(OPTION) cylinder.cpp

rectangle.o: rectangle.h parameter.h ran.h vec.h shape.h rectangle.cpp
	$(CC) -c $(OPTION) rectangle.cpp

ran.o: ran.h ran.cpp
	$(CC) -c $(OPTION) ran.cpp

root6.o: root6.h parameter.h nr.h vec.h root6.cpp
	$(CC) -c $(OPTION) root6.cpp

zrhqr.o: nrutil.h zrhqr.cpp
	$(CC) -c $(OPTION) zrhqr.cpp 

hqr.o: nrutil.h hqr.cpp
	$(CC) -c $(OPTION) hqr.cpp

balanc.o: balanc.cpp
	$(CC) -c $(OPTION) balanc.cpp

nrutil.o: nrutil.cpp
	$(CC) -c $(OPTION) nrutil.cpp

vec.o:  vec.h parameter.h vec.cpp
	$(CC) -c $(OPTION) vec.cpp

parameter.o: parameter.h parameter.cpp
	$(CC) -c $(OPTION) parameter.cpp

clean:
	rm -f *.o
	rm -f *~
	rm -f ellip3d
