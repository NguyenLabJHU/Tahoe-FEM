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
#include "VolumeT.h"
#include "BoxT.h"
#include "OutputSetT.h"
#include "OutPutLatticeT.h"

using namespace Tahoe;

  // Constructor
MeshAtom::MeshAtom(const FCCT& crystal,const BoxT& shape) 
//  MeshAtom::MeshAtom(const CrystalLatticeT& crystal,const VolumeT& shape)
{
  Crystal = new FCCT(crystal);
  Shape = new BoxT(shape);
}


// Create a mesh of atoms. Return ids if i_id = 1, coordinates if icoor = 1 
// connectivities if iconnect = 1 and printout file if iprint = 1.

void MeshAtom::CreateMeshAtom(int i_id,iArrayT* atomid,
			      int icoor,dArray2DT* coords,
			      int iconnect,iArray2DT* connects,			     
			      int iprint,StringT& program_name,
			      StringT& version, StringT& title, 
			      StringT& input_file,
			      IOBaseT::FileTypeT output_format)
{
  
  if (i_id == 1) 
    {
      if( ( (Shape->GetNumberAtoms()) != atomid->Length()))
	throw eSizeMismatch;
      
      atomid = Shape->GetAtomID();
      
    }
  if (icoor == 1) 
    {
      if( (Shape->GetNumberAtoms()) != coords->MajorDim())
	throw eSizeMismatch;
      if( (Shape->GetDimensions()) != coords->MinorDim())
	throw eSizeMismatch;
      
      coords = Shape->GetAtomCoordinates();
    }
  
  if (iconnect == 1) 
    {
      if( (Shape->GetNumberAtoms()) != connects->MajorDim())
	throw eSizeMismatch;
      if( (Shape->GetDimensions()) != connects->MinorDim())
	throw eSizeMismatch;
      
      connects = Shape->GetAtomConnectivities();
    }
  
  
  if (iprint == 1) 
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
  
}
