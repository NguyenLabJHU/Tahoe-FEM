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
   	//create random grain centers
      for (int i=0; i < NumberofGrains; i++) {
         GrainCenters[i]=*(new dArrayT(dim));
         bool doOver;
         do {
            doOver=false;
         
            for (int j=0; j < nSD; j++) {
                //random multiplying factor [0,1)
               double factor = ( (double)rand() / (double)(RAND_MAX + 1));
               GrainCenters[i][j]=factor * (length(j,1)-length(j,0) );
            }//endfor
         	//if the grain centers are too close together, do it again
         	
            for (int k=0; k<i  && !doOver; k++) {
               if ( dArrayT::Distance (GrainCenters[i], GrainCenters[k] ) < (3*lattice_parameter.Max() ) )
                  doOver=true;
            }//endfor
         	
         } while (doOver); //end do
      }//endfor
   	
   }//endconstructor
   
    PolyT::PolyT(int dim, iArrayT cel, 
      dArrayT lattice_parameter,
      iArrayT which_sort, StringT slt,
           iArrayT per, int NumberofGrains) : BoxT(dim, cel, lattice_parameter, which_sort, slt, per) 
   {
   	//set up global variables
      this->NumberofGrains = NumberofGrains;
      GrainCenters.Dimension(NumberofGrains);
      VolType = "POLY";
   	//create random grain centers
      for (int i=0; i < NumberofGrains; i++) {
         GrainCenters[i]=*(new dArrayT(dim));
         bool doOver;
         do {
            doOver=false;
         
            for (int j=0; j < nSD; j++) {
                //random multiplying factor [0,1)
               double factor = ( (double)rand() / (double)(RAND_MAX + 1));
               GrainCenters[i][j]=factor * (length(j,1)-length(j,0) );
            }//endfor
         	//if the grain centers are too close together, do it again
         	
            for (int k=0; k<i  && !doOver; k++) {
               if ( dArrayT::Distance (GrainCenters[i], GrainCenters[k] ) < (3*lattice_parameter.Max() ) )
                  doOver=true;
            }//endfor
         	
         } while (doOver); //end do
      }//endfor
   	
   }//endconstructor



    PolyT::PolyT(const PolyT& source) : BoxT(source) 
   {
     
      NumberofGrains=source.NumberofGrains;
      GrainCenters=source.GrainCenters;
      
   
      VolType = "POLY";
   }


    void PolyT::CreateLattice(CrystalLatticeT* templateLattice) 
   {
      nArrayT<dArrayT> AtomsinGrain[NumberofGrains];
      for (int i=0; i<NumberofGrains; i++) {	AtomsinGrain[i]=*(new nArrayT<dArrayT>);	 }
   			
      for (int currentGrain =0 ; currentGrain < NumberofGrains; currentGrain++) {
      	
      //call superclass CreateLattice here
         BoxT::CreateLattice(GenerateLattice(templateLattice));
         int storedAtoms =0;
         for (int currentAtom = 0; currentAtom < nATOMS; currentAtom++ ){
             //check to see that the atom is closer to the current grain than any other grain
            bool atomIsClosestToCurrentGrain = true;
            dArrayT currentCoord(nSD);
            for (int i=0; i<nSD; i++) {
               currentCoord[i]=atom_coord(currentAtom,i);
            }
         
            for (int tempGrain=0; tempGrain < NumberofGrains &&atomIsClosestToCurrentGrain; tempGrain ++) {
               if ( dArrayT::Distance( currentCoord,GrainCenters[currentGrain]) > dArrayT::Distance (currentCoord, GrainCenters[tempGrain]) ) atomIsClosestToCurrentGrain = false;
            }
               // check to see if atom is overlapping another atom at the boundaries
               //unoptimized
            bool noOverlap=true;
            double tolerance = .35 * templateLattice->GetLatticeParameters().Max() * templateLattice->GetLatticeParameters().Max();
            for (int tempGrain=0; tempGrain<NumberofGrains && noOverlap; tempGrain++) {
               if (tempGrain != currentGrain) {
                  for (int tempAtom=0; tempAtom < AtomsinGrain[currentGrain].Length() && noOverlap; tempAtom ++ ) {
                     if( dArrayT::Distance( currentCoord , AtomsinGrain[currentGrain][tempAtom] ) < tolerance) noOverlap =false;
                  }
               }
            }
          //if everything is good, "keep" the atom
            if (noOverlap && atomIsClosestToCurrentGrain) {
               storedAtoms++;
               AtomsinGrain[currentGrain].Resize(storedAtoms, currentCoord);
            }//endif
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
               if (currentGrain==0) atom_coord.RowCopy(currentAtom, AtomsinGrain[currentGrain][currentAtom]);
               else atom_coord.RowCopy(currentAtom + AtomsinGrain[currentGrain-1].Length(), AtomsinGrain[currentGrain][currentAtom] );
            } //endfor
         }//endfor
      
      //clean up
         delete AtomsinGrain;
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
      
      	
      	   
      }//end for
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
         return new FCCT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),1, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "BCC")
         return new BCCT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),1, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "DIA")
         return new DIAT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),1, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "HEX")
         return new HEXT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),1, rot_matrix, templateLattice->GetAngleRotation() );
      else if(sLATTYPE == "CORUN")
         return new CORUNT(nSD, templateLattice->GetNUCA(), templateLattice->GetLatticeParameters(),1, rot_matrix, templateLattice->GetAngleRotation() );
   		
   }
	
