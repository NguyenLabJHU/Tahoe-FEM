// DEVELOPMENT
#include "MeshAtomT.h"

#include "ExceptionCodes.h"

#include "fstreamT.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iAutoArrayT.h"

#include "CrystalLatticeT.h"
#include "FCCT.h"
#include "BCCT.h"
#include "DIAT.h"
#include "HEXT.h"
#include "CORUNT.h"

#include "VolumeT.h"
#include "BoxT.h"
#include "AsperityT.h"


#include "OutputSetT.h"
#include "OutPutLatticeT.h"
#include "OutputBaseT.h"

   using namespace Tahoe;

// Constructor
    MeshAtomT::MeshAtomT(StringT which_latticetype,int nsd,int nuca,
    	   dArrayT latticeparameter,StringT which_shape,
    	   int whichunit,dArray2DT len, iArrayT cel,int irot,
    	   dArray2DT mat_rot,double angle,iArrayT isort,iArrayT per, int NumberofGrains, int random_seed)
   {
    
     
      if(which_shape == "BOX")
      {
         if(whichunit==1) 
            Shape = new BoxT(nsd,len,latticeparameter,isort,which_latticetype,per);
         else
            Shape = new BoxT(nsd,cel,latticeparameter,isort,which_latticetype,per);
      }
      else if(which_shape == "ASPERITY")
      {
         if(whichunit==1) 
            Shape = new AsperityT(nsd,len,latticeparameter,isort,which_latticetype,per);
         else
            Shape = new AsperityT(nsd,cel,latticeparameter,isort,which_latticetype,per);
      }
      else if (which_shape=="POLY") {
         if (whichunit==1) {
			if (per ==0) len*=2;
			Shape = new PolyT(nsd,len,latticeparameter,isort,which_latticetype,per, NumberofGrains, random_seed);
         }
         else{
			if (per==0)cel*=2;
			Shape = new PolyT(nsd,cel,latticeparameter,isort,which_latticetype,per, NumberofGrains, random_seed);
         }
      }
      else
      {
         cout << "Shape can only be BOX and not: " << which_shape << "\n";
         throw eBadInputValue;
      }
           	 				 
      if(which_latticetype == "FCC")
         Crystal = new FCCT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
      else if(which_latticetype == "BCC")
         Crystal = new BCCT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
      else if(which_latticetype == "DIA")
         Crystal = new DIAT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
      else if(which_latticetype == "HEX")
         Crystal = new HEXT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
      else if(which_latticetype == "CORUN")
         Crystal = new CORUNT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
      else
         throw eBadInputValue;
   
   	
      
      Set = 0;
      IOLattice = 0;
   }
   			



    MeshAtomT:: ~MeshAtomT()
   {
      delete Crystal;
      delete Shape;
      if(Set != 0) delete Set;
      if(IOLattice != 0) delete IOLattice;
   }


    int MeshAtomT::CreateMeshAtom()
   {
      Shape->CreateLattice(Crystal);
      return Shape->GetNumberAtoms();
   }

    double MeshAtomT::Volume_of_Mesh()
   {
      return Shape->GetVolume();
   }

    StringT MeshAtomT::TypeOfVolume()
   {
      return Shape->GetTypeOfVolume();
   }

    dArray2DT MeshAtomT::Length()
   {
      return Shape->GetLength();
   }

    iArrayT MeshAtomT::NumberOFCells()
   {
      return Shape->GetNCells();
   } 

    dArray2DT* MeshAtomT::ReturnCoordinates()
   {   
      return Shape->GetAtomCoordinates();
   }

    dArray2DT* MeshAtomT::ReturnBounds()
   {   
      return Shape->GetAtomBounds();
   }
  
    void MeshAtomT::BuildIOFile(StringT& program_name,
    		   StringT& version, StringT& title, 
    		   StringT& input_file,
    		   IOBaseT::FileTypeT output_format,
    		   iArrayT per)
   {
      Shape->CalculateBounds(Crystal);
   
      if(IOLattice != 0) delete IOLattice;
      IOLattice = new OutPutLatticeT(cout,program_name,version,title,
         	 input_file,output_format);
      ArrayT<StringT> n_labels(1);
      n_labels[0] = "Atom";
      ArrayT<StringT> e_labels(1);
      e_labels[0] = "point";
   
   //if(Set != 0) delete Set;
   //Set=new OutputSetT(GeometryT::kPoint,*(Shape->GetAtomConnectivities()), n_labels);
   
      if(Set != 0) delete Set;
      Set=new OutputSetT(GeometryT::kPoint, *(Shape->GetAtomID()), *(Shape->GetAtomConnect()), 
           n_labels, e_labels, false);
   
      IOLattice->SetCoordinates(*(Shape->GetAtomCoordinates()),(Shape->GetAtomNumber()));
   
      IOLattice->SetBounds(*(Shape->GetAtomBounds()));
      IOLattice->SetTypes(*(Shape->GetAtomTypes()));
      IOLattice->SetParts(*(Shape->GetAtomParts()));
   
   /* Add nodeset ??? (SA Mon Aug  4 2003)
   const iArray2DT nodeset = *((*Shape->GetAtomConnect())[0]);
   const StringT setID = (*Shape->GetAtomID())[0];
   
   const iArrayT nodeset0(Shape->GetNumberOfAtoms());
   for (int p=0;p<nodeset.MajorDim();p++) 
    nodeset0[p] = nodeset(p)[0];
   
   IOLattice->AddNodeSet(nodeset0,setID);
   
   if (Crystal->GetNumberOfType() > 1) 
    {
      const iArray2DT nodeset1 = *((*Shape->GetAtomConnect())[1]);
      const StringT setID1 = (*Shape->GetAtomID())[1];
   
      for (int p=0;p<nodeset1.MajorDim();p++) 
      	nodeset0[p] = nodeset1(p)[0];
   
      IOLattice->AddNodeSet(nodeset0,setID);
    }
   */
   
      IOLattice->AddElementSet(*Set);
   
      IOLattice->WriteGeometry();
   }
