#include <iostream.h>
#include "BoxT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "VolumeT.h"
#include "CrystalLatticeT.h"

#include "FCCT.h"
#include "BCCT.h"
#include "DIAT.h"
#include "HEXT.h"
#include "CORUNT.h"

#include "nArrayT.h"
#include "Rotate3DT.h"
#include "Rotate2DT.h"

   using namespace Tahoe;

    class PolyT : public BoxT 
   {
   	      
   protected:
   
   	//info about the grains in this polycrystalline structure
      int NumberofGrains;
      nArrayT<dArrayT> GrainCenters;
      double MaxGrainSeparation;
      dArray2DT SizeofLattice;
      dArrayT PeriodicLength;
      bool IsPeriodic;
      double Tolerance;
      int kRandomSeed;
      
   
   public:
   
   //Constructor
      PolyT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per, int NumGrains, int random_seed);
   
      PolyT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per, int NumGrains, int random_seed);
   
   //Destructor
       ~PolyT(){};
   
   // Copy constructor
      PolyT(const PolyT& source);
   
      void CreateLattice(CrystalLatticeT* templateLattice); 
   
   protected:
      CrystalLatticeT* GenerateLattice(CrystalLatticeT* templateLattice);
   	dArrayT GenerateRotation();
   	dArrayT GetAtom(Rotate3DT rotationMatrix, dArrayT grainCenter, int AtomNumber);
   	dArrayT GetAtom(Rotate2DT rotationMatrix, dArrayT grainCenter, int AtomNumber);
		void MakeGrains(dArrayT lattice_parameter);
	dArrayT ImposePBC(dArrayT currentCoord, dArrayT GrainCenter);
	bool CheckCurrentGrain(dArrayT currentCoord, int currentGrain, double distanceToCurrentGrain);
	bool CheckNoOverlap (dArrayT currentCoord, nArrayT <nArrayT<dArrayT> > *AtomsinGrain, int currentGrain, CrystalLatticeT *templateLattice);
	void Initialize();
   };


