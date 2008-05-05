RECOMMEND = -O3 -funroll-all-loops -fpeel-loops
OPTIMIZE1 = -fprofile-arcs
OPTIMIZE2 = -fbranch-probalilities

ellip3d: ellip3d.o assembly.o cylinder.o rectangle.o particle.o root6.o ran.o balanc.o zrhqr.o hqr.o vec.o nrutil.o
	g++ $(RECOMMEND) -Wall -gstabs -o ellip3d ellip3d.o assembly.o cylinder.o rectangle.o particle.o root6.o ran.o balanc.o zrhqr.o hqr.o vec.o nrutil.o

ellip3d.o: vec.h gradation.h rectangle.h cylinder.h assembly.h main.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs -o ellip3d.o main.cpp

assembly.o: const.h vec.h gradation.h boundarytgt.h particle.h contact.h boundary.h rectangle.h cylinder.h assembly.h assembly.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs assembly.cpp

particle.o: const.h vec.h gradation.h boundarytgt.h contact.h boundary.h rectangle.h cylinder.h ran.h root6.h particle.h particle.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs particle.cpp

cylinder.o: vec.h shape.h const.h ran.h cylinder.h cylinder.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs cylinder.cpp

rectangle.o: vec.h shape.h ran.h rectangle.h rectangle.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs rectangle.cpp

ran.o: ran.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs ran.cpp

root6.o: const.h vec.h root6.h nr.h root6.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs root6.cpp

zrhqr.o: nrutil.h zrhqr.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs zrhqr.cpp 

hqr.o: nrutil.h hqr.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs hqr.cpp

balanc.o: balanc.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs balanc.cpp

nrutil.o: nrutil.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs nrutil.cpp

vec.o: const.h vec.h vec.cpp
	g++ -c $(RECOMMEND) -Wall -gstabs vec.cpp

clean:
	rm -f *.o
	rm -f *~
	rm -f ellip3d
