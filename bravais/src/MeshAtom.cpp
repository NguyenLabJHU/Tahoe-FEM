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
#include "VolumeT.h"
#include "BoxT.h"
#include "OutputSetT.h"
#include "OutPutLatticeT.h"

using namespace Tahoe;

// Constructor
MeshAtom::MeshAtom(StringT which_latticetype,int nsd,int nuca,
		   dArrayT latticeparameter,StringT which_shape,
		   int whichunit,dArray2DT len, iArrayT cel,
		   dArray2DT mat_rot,double angle)
{
  if(which_latticetype == "FCC")
    Crystal = new FCCT(nsd,nuca,latticeparameter[0],mat_rot,angle);
  else if(which_latticetype == "BCC")
    Crystal = new BCCT(nsd,nuca,latticeparameter[0],mat_rot,angle);
  else if(which_latticetype == "DIA")
    Crystal = new DIAT(nsd,nuca,latticeparameter[0],mat_rot,angle);
  else
    {
      throw eBadInputValue;
    }
  
  if(which_shape == "BOX")
    {
      if(whichunit==1) 
	Shape = new BoxT(nsd,len,latticeparameter);
      else
	Shape = new BoxT(nsd,cel,latticeparameter);
    }
  else
    {
      cout << "Shape can only be BOX and not: " << which_shape << "\n";
      throw eBadInputValue;
    }
  
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

  IOLattice = new OutPutLatticeT(cout,program_name,version,title,
				 input_file,output_format,
				 *(Shape->GetAtomBounds()),
				 *(Shape->GetAtomType()));
  
  ArrayT<StringT> n_labels(1);
  n_labels[0] = "Atom";

  Set=new OutputSetT(*(Shape->GetAtomNames()), GeometryT::kPoint, 
		     *(Shape->GetAtomConnectivities()), n_labels);
  
  IOLattice->SetCoordinates(*(Shape->GetAtomCoordinates()),(Shape->GetAtomID()));
  IOLattice->AddElementSet(*Set);
  
  IOLattice->WriteGeometry();
}

