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
		   int whichunit,dArrayT len_cel,dArrayT rot_vec)
{
  if(which_latticetype == "FCC")
    Crystal = new FCCT(nsd,nuca,latticeparameter[0],rot_vec);
  else if(which_latticetype == "BCC")
    Crystal = new BCCT(nsd,nuca,latticeparameter[0],rot_vec);
  else if(which_latticetype == "DIA")
    Crystal = new DIAT(nsd,nuca,latticeparameter[0],rot_vec);
  else
    {
      throw eBadInputValue;
    }
  
  if(which_shape == "BOX")
    Shape = new BoxT(nsd,whichunit,len_cel,latticeparameter);
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
  Shape->CalculateVolume();
  return Shape->GetVolume();
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
  
void MeshAtom::BuildIOFile(StringT& program_name,
			   StringT& version, StringT& title, 
			   StringT& input_file,
			   IOBaseT::FileTypeT output_format)
{
  IOLattice = new OutPutLatticeT(cout,program_name,version,title,
				 input_file,output_format);
  
  ArrayT<StringT> n_labels(1);
  n_labels[0] = "Atom";
      
  Set=new OutputSetT(*(Shape->GetAtomNames()), GeometryT::kPoint, 
		     *(Shape->GetAtomConnectivities()), n_labels);
  
  IOLattice->SetCoordinates(*(Shape->GetAtomCoordinates()),(Shape->GetAtomID()));
  IOLattice->AddElementSet(*Set);
  
  IOLattice->WriteGeometry();
}
