#include "PolyT.h"
#include "BoxT.h"

#include "ExceptionCodes.h"
#include "iArrayT.h"

#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

#include "FCCT.h"
#include "BCCT.h"
#include "DIAT.h"
#include "HEXT.h"
#include "CORUNT.h"

#include "nArrayT.h"
#include "Rotate3DT.h"
#include "Rotate2DT.h"

#include "time.h"
#include <math.h>
#include <stdlib.h>

    PolyT::PolyT(int dim, dArray2DT len, 
    dArrayT lattice_parameter,
    iArrayT which_sort, StringT slt,
    iArrayT per, int NumberofGrains) : BoxT(dim,len, lattice_parameter, which_sort, slt, per) 
   {
   //set up global variables
      this->NumberofGrains = NumberofGrains;
      GrainCenters.Dimension(NumberofGrains);
      VolType = "POLY";
      SizeofLattice=length;
      if (per == 0) SizeofLattice/=2;
   //create random grain centers
      MaxGrainSeparation=0.0;
      MakeGrains(lattice_parameter);
       
   }//endconstructor
   
    PolyT::PolyT(int dim, iArrayT cel, 
    dArrayT lattice_parameter,
    iArrayT which_sort, StringT slt,
    iArrayT per, int NumberofGrains) : BoxT(dim, cel, lattice_parameter, which_sort, slt, per) 
   {
   //set up global variables
      this->NumberofGrains = NumberofGrains;
      GrainCenters.Dimension(NumberofGrains);
      SizeofLattice=length;
      if(per ==0) SizeofLattice/=2;
      VolType = "POLY";
      MaxGrainSeparation=0.0;
      MakeGrains(lattice_parameter);
      cout<<"nATOMS= " << nATOMS <<"\n";
      cout<<"SizeofLattice: "<<SizeofLattice<<"\n";
   
   }//endconstructor



    PolyT::PolyT(const PolyT& source) : BoxT(source) 
   {
     
      NumberofGrains=source.NumberofGrains;
      GrainCenters=source.GrainCenters;
      MaxGrainSeparation=source.MaxGrainSeparation;
      VolType = "POLY";
      SizeofLattice = source.SizeofLattice;
   }


    void PolyT::CreateLattice(CrystalLatticeT* templateLattice) 
   {
   
//      dArray2DT virtualGrainBounds(nSD,2);
//   
//      for (int i = 0; i < nSD; i++) {
//         virtualGrainBounds(i, 0)= -MaxGrainSeparation/2.0;
//         virtualGrainBounds(i, 1)= MaxGrainSeparation/2.0;
//      }
//     //make a virtual grain
//      BoxT::BoxT(nSD, virtualGrainBounds, templateLattice->GetLatticeParameters(), WhichSort, sLATTYPE, pbc);

      BoxT::CreateLattice(templateLattice);
     
     //initialize an array to store the atom coordinates associated with each grain   
      nArrayT <nArrayT<dArrayT> > AtomsinGrain(NumberofGrains);
   
      nArrayT <nArrayT<int> > TypesinGrain(NumberofGrains);  
      for (int i=0; i<NumberofGrains; i++) {
         AtomsinGrain[i]=*(new nArrayT<dArrayT>);
         TypesinGrain[i]=*(new nArrayT<int>);
      }
         
      for (int currentGrain =0 ; currentGrain < NumberofGrains; currentGrain++) {
         dArrayT rotation_vector=GenerateRotation();
         int storedAtoms =0;
         Rotate3DT rotate3;
         Rotate2DT rotate2;
         if (nSD==2) rotate2 = *(new Rotate2DT(rotation_vector[0])); 
         if (nSD==3) rotate3= *(new Rotate3DT(rotation_vector[0], rotation_vector[1], rotation_vector[2]));
         for (int currentAtom = 0; currentAtom < nATOMS; currentAtom++ ){
            dArrayT currentCoord(nSD);
			//grab the coordinates of the current atom;
            if (nSD==2) currentCoord = GetAtom(rotate2, GrainCenters[currentGrain], currentAtom );
            if (nSD==3) currentCoord = GetAtom(rotate3, GrainCenters[currentGrain], currentAtom );

			//get the distance to the current grain center and impose PeriodicBCs
			double distanceToCurrentGrain = dArrayT::Distance(currentCoord, GrainCenters[currentGrain]);
			if(pbc!=0) currentCoord=ImposePBC(currentCoord, GrainCenters[currentGrain]);
			//check to see if atom is outside the bounds of the lattice region
            bool atomIsInside = true;
            for (int tempDim=0; tempDim < nSD && atomIsInside; tempDim++) {
               atomIsInside = currentCoord[tempDim] >= SizeofLattice(tempDim,0) &&currentCoord[tempDim] <= SizeofLattice(tempDim,1);
            }        
            
         	//check to see that the atom is closer to the current grain than any other grain
         	bool atomIsClosestToCurrentGrain=false;
			if (atomIsInside) atomIsClosestToCurrentGrain = CheckCurrentGrain(currentCoord, currentGrain, distanceToCurrentGrain);
            

         
         	// check to see if atom is overlapping another atom at the boundaries
            bool noOverlap=false;
            if(atomIsInside&&atomIsClosestToCurrentGrain) noOverlap=CheckNoOverlap(currentCoord, &AtomsinGrain, currentGrain, templateLattice);

         //if everything is good, "keep" the atom
            if (noOverlap && atomIsClosestToCurrentGrain && atomIsInside) {
               storedAtoms++;
            
               AtomsinGrain[currentGrain].Resize(storedAtoms, currentCoord);
               TypesinGrain[currentGrain].Resize(storedAtoms, atom_types[currentAtom]);
            
            }//endif
         }//endfor
         cout <<"stored atoms: " <<storedAtoms<<"\n";
      }//endfor
         //now that we've gone through all the grains, copy these temporary atoms back into atom_coords
         //first, calculate total number of atoms
      nATOMS=0;
      for (int currentGrain = 0; currentGrain < NumberofGrains; currentGrain++) 
         nATOMS += AtomsinGrain[currentGrain].Length();
      cout<<"nATOMS ="<<nATOMS <<"\n";
      atom_coord.Free();
      atom_coord.Dimension(nATOMS, nSD);
      int totalAtoms=0;
      for (int currentGrain = 0; currentGrain<NumberofGrains; currentGrain++) {
         for (int currentAtom=0; currentAtom< AtomsinGrain[currentGrain].Length();currentAtom++) {
            atom_coord.SetRow(totalAtoms, AtomsinGrain[currentGrain][currentAtom]);
            atom_types[totalAtoms]=TypesinGrain[currentGrain][currentAtom];
            totalAtoms++;
         
         } //endfor
      }//endfor
      cout<<"finished copying atoms\n";
      //clean up
      AtomsinGrain.Free();
   
   // Get atoms coordinates
      atom_number.Dimension(nATOMS);
      
   
      for(int m=0; m < nATOMS ; m++) 
      {
         atom_number[m] = m;
      
      }
   
      atom_names = "Poly";
   
      // Update lengths
      length = SizeofLattice;
      cout<<"length updated \n";
      //Calculate volume here
      switch(nSD) {
         case 2:
            volume = (length(0,1)-length(0,0))*(length(1,1)-length(1,0));
            break;
         case 3:
            volume = (length(0,1)-length(0,0))*         
                  (length(1,1)-length(1,0))*
                  (length(2,1)-length(2,0));
            break;
      }//endswitch
      
      // Sort Lattice 
      //if(WhichSort  != 0) SortLattice(templateLattice);
      int ntype = templateLattice->GetNTYPE();
      // Create types and connectivities
      if(ntype > 2) 
         cout << "WARNING: ** nTypes == 2  maximum for ensight output ** \n";
      
      atom_ID.Dimension(ntype);
      atom_connect.Dimension(ntype);
      
      type1.Dimension(nATOMS,1);type1 = 0;
      type2.Dimension(nATOMS,1);type2 = 0;
      
      int n=0,m=0;
      for (int p=0;p<nATOMS;p++)
      {
         if ( atom_types[p] == 1) 
         {
            type1(n)[0] = p;
            n++;
         }	
         else
         {
            type2(m)[0] = p;
            m++;
         }	
      }//endfor
      
      atom_ID[0] = "type1";
      if (ntype > 1) atom_ID[1] = "type2";
      
      if (n < nATOMS && n > 0) type1.Resize(n);
      if (m < nATOMS && m > 0) type2.Resize(m);
      atom_connect[0] = &type1;
      if (ntype > 1) atom_connect[1] = &type2;
      
      
      // Create parts
      atom_parts.Dimension(nATOMS);
      atom_parts = 1;
      cout<<"finished lattice creation \n";
    
   }//end function



    dArrayT PolyT::GenerateRotation() {
   
      dArrayT rot_vector (nSD);
      
   
      if (nSD==3) {
      //create a random rotation matrix
      //random euler angles
         double phi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX));
         double theta = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX));
         double psi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX));
         rot_vector[0]=phi;
         rot_vector[1]=theta;
         rot_vector[2]=psi;
      }
      if (nSD==2) {
         double angle=3.14159265/2.0 * ( (double)rand() / (double)(RAND_MAX ));
         rot_vector[0]=angle;
      }
      return rot_vector;
   }

    dArrayT PolyT::GetAtom(Rotate3DT rotationMatrix, dArrayT grainCenter, int AtomNumber) {
      dArrayT atomCoord(nSD);
      for (int i=0;i<nSD;i++) {
         atomCoord[i]=atom_coord(AtomNumber, i);
      }
      dArrayT newCoord(nSD);
      newCoord=rotationMatrix.RotateVectorIn(atomCoord);
   
      newCoord+=grainCenter;
      return newCoord;
   }
	
    dArrayT PolyT::GetAtom(Rotate2DT rotationMatrix, dArrayT grainCenter, int AtomNumber) {
      dArrayT atomCoord(nSD);
      for (int i=0;i<nSD;i++) {
         atomCoord[i]=atom_coord(AtomNumber, i);
      }
      atomCoord=rotationMatrix.RotateVectorIn(atomCoord);
      atomCoord+=grainCenter;
      return atomCoord;
   }

    void PolyT::MakeGrains (dArrayT lattice_parameter) {
    	srand(time(0));
   	
      for (int i=0; i < NumberofGrains; i++) {

         GrainCenters[i]=*(new dArrayT(nSD));
         bool doOver;
         do {
            doOver=false;
         
            for (int j=0; j < nSD; j++) {
            //random multiplying factor [0,1) to randomly place the grain center
               double factor = ( (double)rand() / (double)(RAND_MAX));
               cout <<"factor: "<<factor<<"\n";
               GrainCenters[i][j]=factor * (SizeofLattice(j,1)-SizeofLattice(j,0) ) + SizeofLattice(j,0);
            }//endfor
         
         //if the grain centers are too close together, do it again
         
            for (int k=0; k<i  && !doOver; k++) {
               if ( dArrayT::Distance (GrainCenters[i], GrainCenters[k] ) < (3*lattice_parameter.Max() ) )
                  doOver=true;
            }//endfor
         	
         } while (doOver); //end do
         cout<<"Grain "<<i<<": "<<GrainCenters[i]<<"\n";
      //calculate the maximum grain separation
         for (int tempGrain=0; tempGrain<i; tempGrain++) {
            if (dArrayT::Distance(GrainCenters[i],GrainCenters[tempGrain]) > MaxGrainSeparation ) 
               MaxGrainSeparation = dArrayT::Distance(GrainCenters[i],GrainCenters[tempGrain]);
         	
         }
         for(int dim=0;dim<nSD;dim++) {
            if (GrainCenters[i][dim]-SizeofLattice(dim,0) > MaxGrainSeparation ) MaxGrainSeparation = GrainCenters[i][dim]-SizeofLattice(dim,0);
            if (SizeofLattice(dim,1)-GrainCenters[i][dim] > MaxGrainSeparation ) MaxGrainSeparation = SizeofLattice(dim,1)-GrainCenters[i][dim];
         }
      }//endfor
      cout<<"Max Grain Separation: "<<MaxGrainSeparation<<"\n";
   }
   
dArrayT PolyT::ImposePBC(dArrayT currentCoord, dArrayT GrainCenter){
	for (int dim=0; dim<nSD ; dim++){
		if (currentCoord[dim]<SizeofLattice(dim,0) && pbc[dim]==1) currentCoord[dim]+= SizeofLattice(dim,1)*2.0;
		else if (currentCoord[dim]>SizeofLattice(dim,1) && pbc[dim]==1) currentCoord[dim]-=SizeofLattice(dim,1)*2.0;
	}
	return currentCoord;
}

bool PolyT::CheckCurrentGrain(dArrayT currentCoord, int currentGrain, double distanceToCurrentGrain){
	bool atomIsCloserToCurrentGrain=true;
	for (int tempGrain=0; tempGrain < NumberofGrains && atomIsCloserToCurrentGrain; tempGrain++) {
		if (tempGrain != currentGrain) {
			dArrayT distance(nSD);
			distance.DiffOf(currentCoord , GrainCenters[tempGrain]);
			for (int dim=0; dim<nSD &&pbc!=0; dim++) {
				if (distance[dim]< SizeofLattice(dim,0) && pbc[dim]==1 )distance[dim]+= SizeofLattice(dim,1)*2.0;
				else if (distance[dim] > SizeofLattice(dim,1) && pbc[dim]==1) distance[dim]-=SizeofLattice(dim,1)*2.0;
			}
			atomIsCloserToCurrentGrain = distanceToCurrentGrain < distance.Magnitude();
		}
	}
	return atomIsCloserToCurrentGrain;
}

bool PolyT::CheckNoOverlap (dArrayT currentCoord, nArrayT <nArrayT<dArrayT> > *AtomsinGrain, int currentGrain, CrystalLatticeT *templateLattice) {
    double tolerance = .35 * templateLattice->GetLatticeParameters().Max() * templateLattice->GetLatticeParameters().Max();
    bool noOverlap=true;
    for (int tempGrain=0; tempGrain<currentGrain && noOverlap; tempGrain++) {
		for (int tempAtom=0; tempAtom < (*AtomsinGrain)[tempGrain].Length() && noOverlap; tempAtom++ ) {
			dArrayT distance(nSD);
			distance.DiffOf(currentCoord,(*AtomsinGrain)[tempGrain][tempAtom]);
			for (int dim=0; dim<nSD && pbc!=0; dim++) {
				if (distance[dim]< SizeofLattice(dim,0) && pbc[dim]==1 )distance[dim]+= SizeofLattice(dim,1)*2.0;
				else if (distance[dim] > SizeofLattice(dim,1) && pbc[dim]==1) distance[dim]-=SizeofLattice(dim,1)*2.0;
			}
			noOverlap = distance.Magnitude() >= tolerance;
       }
    }
    return noOverlap;
}