#ifndef _ATOMMESH_H_
#define _ATOMMESH_H_

#include <iostream>
#include <fstream.h>

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iAutoArrayT.h"
#include "IOBaseT.h"
#include "GeometryT.h"

#include "ifstreamT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

#include "CrystalLatticeT.h"
#include "FCCT.h"
#include "VolumeT.h"
#include "BoxT.h"
#include "OutputSetT.h"
#include "OutPutLatticeT.h"


namespace Tahoe {

/** Class to construct a mesh of atoms defined by
 *  its crytallography and its shape. 
 *  Create an output file in different formats.
 *  
 *  Note: To create a mesh with another types than
 *        FCC or Box, declare another constructor
 *        with the right parameters. 
 *        Example: to create an FCC lattice in a box, 
 *                 call   MeshAtom(const FCCT& crystal,
 *  	                           const BoxT& shape,
 *	                           const OutputSetT& outset); 
 *				   
 *                 to create a BCC  lattice in a sphere,
 *                 call   MeshAtom(const BCCT& crystal,
 *  	                           const SphereT& shape,
 *	                           const OutputSetT& outset); 
 *
 **/

class MeshAtom {

 protected:

  CrystalLatticeT* Crystal;
  VolumeT* Shape;
  OutputSetT* Set;
  OutPutLatticeT* IOLattice;

 public:

  // Constructor
    MeshAtom(const FCCT& crystal,const BoxT& shape); 
    //MeshAtom(const CrystalLatticeT& crystal,const VolumeT& shape); 
  
  // Destructor
  ~MeshAtom(){};

// Create a mesh of atoms. Return ids if i_id = 1, coordinates if icoor = 1 
// connectivities if iconnect = 1 and printout file if iprint = 1.
void MeshAtom::CreateMeshAtom(int i_id,iArrayT* atomid,
			      int icoor,dArray2DT* coords,
			      int iconnect,iArray2DT* connects,			     
			      int iprint,StringT& program_name,
			      StringT& version, StringT& title, 
			      StringT& input_file,
			      IOBaseT::FileTypeT output_format);

};

} // namespace Tahoe 


#endif

    
