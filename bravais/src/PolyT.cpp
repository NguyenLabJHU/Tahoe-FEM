#include "PolyT.h"
#include "BoxT.h"

#include "ExceptionCodes.h"
#include "ifstreamT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
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
#include <math.h>
#include <stdlib.h>

    PolyT::PolyT(int dim, dArray2DT len, 
    dArrayT lattice_parameter,
    iArrayT which_sort, StringT slt,
    iArrayT per, int NumberofGrains) : BoxT(dim, len, lattice_parameter, which_sort, slt, per) 
   {
   //set up global variables
      this->NumberofGrains = NumberofGrains;
      GrainCenters.Dimension(NumberofGrains);
      VolType = "POLY";
      SizeofLattice=length;
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
      VolType = "POLY";
      MaxGrainSeparation=0.0;
      MakeGrains(lattice_parameter);
   
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
      dArray2DT virtualGrainBounds(nSD,2);
      for (int i = 0; i < nSD; i++) {
         virtualGrainBounds(i, 0)= -MaxGrainSeparation/2.0;
         virtualGrainBounds(i, 1)= MaxGrainSeparation/2.0;
      }
     //make a virtual grain
      BoxT::BoxT(nSD, virtualGrainBounds, templateLattice->GetLatticeParameters(), WhichSort, sLATTYPE, pbc);
      BoxT::CreateLattice(templateLattice);
     
     //initialize an array to store the atom coordinates associated with each grain   
      nArrayT<dArrayT> AtomsinGrain[NumberofGrains];
      nArrayT<int> TypesinGrain[NumberofGrains];  
      for (int i=0; i<NumberofGrains; i++) {
         AtomsinGrain[i]=*(new nArrayT<dArrayT>);
         TypesinGrain[i]=*(new nArrayT<int>);
      }
         
      for (int currentGrain =0 ; currentGrain < NumberofGrains; currentGrain++) {
         dArrayT rotation_vector=GenerateRotation();
         int storedAtoms =0;
         Rotate3DT rotate3(rotation_vector[0], rotation_vector[1], rotation_vector[2]);
         Rotate2DT rotate2(rotation_vector[0]);
         for (int currentAtom = 0; currentAtom < nATOMS; currentAtom++ ){
            dArrayT currentCoord(nSD);
         //grab the coordinates of the current atom;
            if (nSD==2) currentCoord = GetAtom(rotate2, GrainCenters[currentGrain], currentAtom );
            if (nSD==3) currentCoord = GetAtom(rotate3, GrainCenters[currentGrain], currentAtom );
         //check to see that the atom is closer to the current grain than any other grain
            bool atomIsClosestToCurrentGrain = true;
            for (int tempGrain=0; tempGrain < NumberofGrains &&atomIsClosestToCurrentGrain; tempGrain ++) {
               if ( dArrayT::Distance( currentCoord,GrainCenters[currentGrain]) > dArrayT::Distance (currentCoord, GrainCenters[tempGrain]) ) atomIsClosestToCurrentGrain = false;
            }
         //check to see if atom is outside the bounds of the lattice region
            bool atomIsInside = true;
            for (int tempDim=0; tempDim < nSD && atomIsInside; tempDim++) {
               atomIsInside = currentCoord[tempDim] >= SizeofLattice(tempDim,0) &&currentCoord[tempDim] <= SizeofLattice(tempDim,1);
            }
         
         // check to see if atom is overlapping another atom at the boundaries
            bool noOverlap=true;
            double tolerance = .35 * templateLattice->GetLatticeParameters().Max() * templateLattice->GetLatticeParameters().Max();
            for (int tempGrain=0; tempGrain<currentGrain && noOverlap; tempGrain++) {
               for (int tempAtom=0; tempAtom < AtomsinGrain[tempGrain].Length() && noOverlap; tempAtom ++ ) {
                  if( dArrayT::Distance( currentCoord , AtomsinGrain[tempGrain][tempAtom] ) < tolerance) noOverlap =false;
               }
            }
         //if everything is good, "keep" the atom
            if (noOverlap && atomIsClosestToCurrentGrain && atomIsInside) {
               storedAtoms++;
               AtomsinGrain[currentGrain].Resize(storedAtoms, currentCoord);
               TypesinGrain[currentGrain].Resize(storedAtoms, atom_types[currentAtom]);
            }//endif
         }//endfor
         
      }//endfor
         //now that we've gone through all the grains, copy these temporary atoms back into atom_coords
         //first, calculate total number of atoms
      nATOMS=0;
      for (int currentGrain = 0; currentGrain < NumberofGrains; currentGrain++) 
         nATOMS += AtomsinGrain[currentGrain].Length();
      atom_coord.Free();
      atom_coord.Dimension(nATOMS, nSD);
      
      for (int currentGrain = 0; currentGrain<NumberofGrains; currentGrain++) {
         for (int currentAtom=0; currentAtom< AtomsinGrain[currentGrain].Length();currentAtom++) {
            if (currentGrain==0) {
               atom_coord.RowCopy(currentAtom, AtomsinGrain[currentGrain][currentAtom]);
               atom_types[currentAtom]=TypesinGrain[currentGrain][currentAtom];
            }
            else {
               atom_coord.RowCopy(currentAtom + AtomsinGrain[currentGrain-1].Length(), AtomsinGrain[currentGrain][currentAtom] );
               atom_types[currentAtom+AtomsinGrain[currentGrain-1].Length()]= TypesinGrain[currentGrain][currentAtom];
            }
         } //endfor
      }//endfor
      
      //clean up
      delete AtomsinGrain;
      
   
   // Get atoms coordinates
      atom_number.Dimension(nATOMS);
      
   
      for(int m=0; m < nATOMS ; m++) 
      {
         atom_number[m] = m;
      }
   
      atom_names = "Poly";
   
      // Update lengths
      length = ComputeMinMax();
      
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
      if(WhichSort  != 0) SortLattice(templateLattice);
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
    
   }//end function


    CrystalLatticeT* PolyT:: GenerateLattice (CrystalLatticeT* templateLattice) {
      dArray2DT rot_matrix = *(new dArray2DT(nSD, nSD));
   	//create a random rotation matrix
   	//random euler angles
      double phi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
      double theta = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
      double psi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
      //define the rotation matrix
      rot_matrix(0,0)=cos(psi)*sin(phi) - cos(theta)*sin(phi)*sin(psi);
      rot_matrix(0,1) = cos(psi)*sin(phi)+ cos(theta)*cos(phi)*sin(psi);
      rot_matrix(0,2)=sin(theta)*sin(psi);
      rot_matrix(1,0)=-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
      rot_matrix(1,1)= -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
      rot_matrix(1,2) = sin(theta)*cos(psi);
      rot_matrix(2,0) = sin(theta)*sin(phi);
      rot_matrix(2,1)= -sin(theta)*cos(phi);
      rot_matrix(2,2) =cos(theta);
        	 				 
      if(sLATTYPE == "FCC")
         return new FCCT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),0, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "BCC")
         return new BCCT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),0, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "DIA")
         return new DIAT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),0, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "HEX")
         return new HEXT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),0, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "CORUN")
         return new CORUNT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),0, rot_matrix, templateLattice->GetAngleRotation() );
   		
   }

    dArrayT PolyT::GenerateRotation() {
      dArrayT rot_vector (nSD);
      if (nSD==3) {
      //create a random rotation matrix
      //random euler angles
         double phi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
         double theta = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
         double psi = 3.14159265 / 2.0 * ( (double)rand() / (double)(RAND_MAX + 1));
         rot_vector[0]=phi;
         rot_vector[1]=theta;
         rot_vector[2]=psi;
      }
      if (nSD==2) {
         double angle=3.14159265 * ( (double)rand() / (double)(RAND_MAX + 1));
         rot_vector[0]=angle;
      }
      return rot_vector;
   }

    dArrayT PolyT::GetAtom(Rotate3DT rotationMatrix, dArrayT grainCenter, int AtomNumber) {
      dArrayT atomCoord(nSD);
      for (int i=0;i<nSD;i++) {
         atomCoord[i]=atom_coord(AtomNumber, i);
      }
      atomCoord=rotationMatrix.RotateVectorIn(atomCoord);
      atomCoord+=grainCenter;
      return atomCoord;
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
      for (int i=0; i < NumberofGrains; i++) {
         GrainCenters[i]=*(new dArrayT(nSD));
         bool doOver;
         do {
            doOver=false;
         
            for (int j=0; j < nSD; j++) {
            //random multiplying factor [0,1) to randomly place the grain center
               double factor = ( (double)rand() / (double)(RAND_MAX + 1));
               GrainCenters[i][j]=factor * (length(j,1)-length(j,0) );
            }//endfor
         
         //if the grain centers are too close together, do it again
         
            for (int k=0; k<i  && !doOver; k++) {
               if ( dArrayT::Distance (GrainCenters[i], GrainCenters[k] ) < (3*lattice_parameter.Max() ) )
                  doOver=true;
            }//endfor
         	
         } while (doOver); //end do
      //calculate the maximum grain separation
         for (int tempGrain=0; tempGrain<i; tempGrain++) {
            if (dArrayT::Distance(GrainCenters[i],GrainCenters[tempGrain]) > MaxGrainSeparation ) 
               MaxGrainSeparation = dArrayT::Distance(GrainCenters[i],GrainCenters[tempGrain]);
         	
         }
      }//endfor
   }