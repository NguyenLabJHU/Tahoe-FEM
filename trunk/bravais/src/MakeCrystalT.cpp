/* $Id: MakeCrystalT.cpp,v 1.3 2002-07-24 01:14:59 saubry Exp $ */


/* Build a mesh of atoms with a format-independent output
 * (saubry, Tue Jul 23 15:31 UTC 2002)
*/


#include "MakeCrystalT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "VolumeT.h"
#include "BoxT.h"
#include "CrystalLatticeT.h"
#include "FCCT.h"

#include "OutputSetT.h"
#include "OutputBaseT.h"
#include "OutPutLatticeT.h"

#include "MeshAtom.h"

void MakeCrystalT::Run() {

  // Read and store data from "data" file	
  StringT inputfile;
  inputfile = "data";
  
  ifstreamT in('%');
  in.open(inputfile);

  int nsd; 
  in >> nsd;
  cout << "Dimension:" << nsd << "\n";

  int whichunit;
  dArrayT length;
  length.Dimension(nsd);

  in >> whichunit;
  if(nsd==2) 
    in >> length[0] >>length[1];
  else if (nsd==3)
    in >> length[0] >>length[1] >> length[2];

  StringT latticetype;
  in >> latticetype;
  cout << "Lattice Type:" << latticetype <<  "\n";
  
  int a,b;
  if (latticetype=="FCC") {
    if(nsd ==2) {a=2;b=2;}
    if(nsd ==3) {a=3;b=4;}
  }
  else {
    cout << "Lattice type has to be FCC...\n";
    throw eBadInputValue;
  }

  double alat;
  in >> alat;

  // Define a crystal lattice: class CrystalLatticeT, FCCT	
  FCCT* crystal_lattice;
  crystal_lattice = new FCCT(a,b,alat);
	
  cout << "Lattice parameter\n" << crystal_lattice->GetLatticeParameters() 
       << "\n";
  crystal_lattice->CalculateDensity();

  // Define a shape: class VolumeT, BoxT
  StringT shape;
  in >> shape;
  cout << "Shape of the domain:" << shape <<  "\n";

  BoxT* crystal_shape;
  if (shape=="box")
    crystal_shape=new BoxT(nsd,whichunit,length,
			crystal_lattice->GetLatticeParameters());
  else 
    {
      cout << "only shape = box is defined so far..\n";
      throw eBadInputValue;
    }

  switch(crystal_shape->GetDimensions()) {
  case 2:
    cout << "Lengths of box: " 
	 << (crystal_shape->GetLength())[0] << " " 
	 << (crystal_shape->GetLength())[1] << "\n";
    break;
  case 3:
    cout << "Lengths of box: " 
	 << (crystal_shape->GetLength())[0] << " " 
	 << (crystal_shape->GetLength())[1] << " " 
	 << (crystal_shape->GetLength())[2] << "\n";
    break;
  default :
    throw eBadInputValue;
  }

  crystal_shape->CalculateVolume();
  cout << "Volume equals " 
       << crystal_shape->GetVolume() << " cubic angstroms" << endl;
  
  crystal_shape->CreateLattice(crystal_lattice);
  cout << "Number of Atoms:" << crystal_shape->GetNumberAtoms() <<  "\n";
 
  in.close();
  
  //Define Mesh
  MeshAtom mesh_atom(*crystal_lattice,*crystal_shape);

  StringT program = "bravais";
  StringT version = "v1.0";
  StringT title = "Lattice for Atoms";
  StringT input = "example";
  
  mesh_atom.CreateMeshAtom(0,NULL,0,NULL,0,NULL,
			   1,program,version,title,input,
			   IOBaseT::kEnSight);
  

}


