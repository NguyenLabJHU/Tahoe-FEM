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
		   dArray2DT mat_rot,double angle,iArrayT isort,iArrayT per)
{
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
    {
      throw eBadInputValue;
    }

  StringT slt = Crystal->GetSLATTYPE();
  
  if(which_shape == "BOX")
    {
      if(whichunit==1) 
	Shape = new BoxT(nsd,len,latticeparameter,isort,slt,per);
      else
	Shape = new BoxT(nsd,cel,latticeparameter,isort,slt,per);
    }
  else if(which_shape == "ASPERITY")
    {
      if(whichunit==1) 
	Shape = new AsperityT(nsd,len,latticeparameter,isort,per);
      else
	Shape = new AsperityT(nsd,cel,latticeparameter,isort,per);
    }
  else
    {
      cout << "Shape can only be BOX and not: " << which_shape << "\n";
      throw eBadInputValue;
    }
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

dArray2DT MeshAtomT::Length()
{
  return Shape->GetLength();
}

iArrayT MeshAtomT::NumberOFCells()
{
  return Shape->GetNCells();
}


iArrayT* MeshAtomT::ReturnAtomID()
{
  return Shape->GetAtomID();
}  

dArray2DT* MeshAtomT::ReturnCoordinates()
{   
  return Shape->GetAtomCoordinates();
}
  
iArray2DT* MeshAtomT::ReturnConnectivities()
{
  return Shape->GetAtomConnectivities();
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

  if(Set != 0) delete Set;
  Set=new OutputSetT(GeometryT::kPoint,*(Shape->GetAtomConnectivities()), n_labels);

  IOLattice->SetCoordinates(*(Shape->GetAtomCoordinates()),(Shape->GetAtomID()));

  IOLattice->SetBounds(*(Shape->GetAtomBounds()));
  IOLattice->SetTypes(*(Shape->GetAtomTypes()));
  IOLattice->SetParts(*(Shape->GetAtomParts()));
 
  IOLattice->AddElementSet(*Set);
  
  IOLattice->WriteGeometry();
}

