// DEVELOPMENT
#include "MeshAtom.h"

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
#include "VolumeT.h"
#include "BoxT.h"
#include "OutputSetT.h"
#include "OutPutLatticeT.h"

using namespace Tahoe;

// Constructor
MeshAtom::MeshAtom(StringT which_latticetype,int nsd,int nuca,
		   dArrayT latticeparameter,StringT which_shape,
		   int whichunit,dArray2DT len, iArrayT cel,int irot,
		   dArray2DT mat_rot,double angle,iArrayT isort)
{
  if(which_latticetype == "FCC")
    Crystal = new FCCT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
  else if(which_latticetype == "BCC")
    Crystal = new BCCT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
  else if(which_latticetype == "DIA")
    Crystal = new DIAT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
  else if(which_latticetype == "HEX")
    Crystal = new HEXT(nsd,nuca,latticeparameter,irot,mat_rot,angle);
  else
    {
      throw eBadInputValue;
    }
  
  if(which_shape == "BOX")
    {
      if(whichunit==1) 
	Shape = new BoxT(nsd,len,latticeparameter,isort);
      else
	Shape = new BoxT(nsd,cel,latticeparameter,isort);
    }
  else
    {
      cout << "Shape can only be BOX and not: " << which_shape << "\n";
      throw eBadInputValue;
    }
  Set = 0;
  IOLattice = 0;
}

MeshAtom:: ~MeshAtom()
{
  delete Crystal;
  delete Shape;
  if(Set != 0) delete Set;
  if(IOLattice != 0) delete IOLattice;
}


int MeshAtom::CreateMeshAtom()
{
  Shape->CreateLattice(Crystal);
  return Shape->GetNumberAtoms();
}

double MeshAtom::Volume_of_Mesh()
{
  return Shape->GetVolume();
}

dArray2DT MeshAtom::Length()
{
  return Shape->GetLength();
}


iArrayT* MeshAtom::ReturnAtomID()
{
  return Shape->GetAtomID();
}  

dArray2DT* MeshAtom::ReturnCoordinates()
{   
  return Shape->GetAtomCoordinates();
}
  
iArray2DT* MeshAtom::ReturnConnectivities()
{
  return Shape->GetAtomConnectivities();
}

dArray2DT* MeshAtom::ReturnBounds()
{   
  return Shape->GetAtomBounds();
}
  
void MeshAtom::BuildIOFile(StringT& program_name,
			   StringT& version, StringT& title, 
			   StringT& input_file,
			   IOBaseT::FileTypeT output_format,
			   iArrayT per)
{
  Shape->CalculateBounds(per,Crystal);
  Shape->CalculateType();

  if(IOLattice != 0) delete IOLattice;

  IOLattice = new OutPutLatticeT(cout,program_name,version,title,
				 input_file,output_format,
				 *(Shape->GetAtomBounds()),
				 *(Shape->GetAtomType()));
  ArrayT<StringT> n_labels(1);
  n_labels[0] = "Atom";

  if(Set != 0) delete Set;

  Set=new OutputSetT(GeometryT::kPoint, 
		     *(Shape->GetAtomConnectivities()), n_labels);
  
  IOLattice->SetCoordinates(*(Shape->GetAtomCoordinates()),(Shape->GetAtomID()));
  IOLattice->AddElementSet(*Set);
  
  cout << "Writing geometry...\n";
  IOLattice->WriteGeometry();
}

