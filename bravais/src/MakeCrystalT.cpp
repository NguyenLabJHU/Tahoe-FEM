/* $Id: MakeCrystalT.cpp,v 1.4 2002-07-24 23:14:56 saubry Exp $ */


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
#include "BCCT.h"
#include "DIAT.h"

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
  dArrayT len_cell;
  len_cell.Dimension(nsd);

  in >> whichunit;
  if(nsd==2) 
    in >> len_cell[0] >>len_cell[1];
  else if (nsd==3)
    in >> len_cell[0] >>len_cell[1] >> len_cell[2];

  StringT latticetype;
  in >> latticetype;
  cout << "Lattice Type:" << latticetype <<  "\n";
  
  int b;
  double alat;
  in >> alat;

  if (latticetype=="FCC") 
    {
      if(nsd == 2) b=2;
      if(nsd == 3) b=4;
    }
  else if (latticetype=="BCC") 
    {
      if(nsd == 2) b=1;
      if(nsd == 3) b=2;
    }
  else if (latticetype=="DIA") 
    {
      if(nsd == 2) b=2;
      if(nsd == 3) b=8;
    }
  else 
    {
      cout << "Lattice type has to be FCC or BCC...\n";
      throw eBadInputValue;
    }
	
  cout << "Lattice parameter: " << alat << "\n";

  StringT shape;
  in >> shape;
  cout << "Shape of the domain:" << shape <<  "\n";

  in.close();
  
  //Define Mesh

  dArrayT latticeparameter;
  latticeparameter.Dimension(nsd);

  for (int i=0; i<nsd; i++)   
    latticeparameter[i] = alat;

  MeshAtom mesh_atom(latticetype,nsd,b,latticeparameter,
		     shape,whichunit,len_cell);

  StringT program = "bravais";
  StringT version = "v1.0";
  StringT title = "Lattice for Atoms";
  StringT input = "example";

  int nb_atoms;
  cout << "\nCreating mesh of atom...\n";
  nb_atoms = mesh_atom.CreateMeshAtom();
  cout << nb_atoms << " atoms in mesh\n";


  cout << "The total volume of the domain is " 
       << mesh_atom.Volume_of_Mesh()
       << " in Angstroms^3 \n";

  dArray2DT coords;
  coords.Dimension(nb_atoms,nsd);

  coords = *(mesh_atom.ReturnCoordinates());
  cout << "Coordinates:\n";
  for (int j=0; j<nb_atoms; j++) 
    cout << coords(j,0) <<  "  " 
	 << coords(j,1) <<  "  " 
	 << coords(j,2) << "\n";
  


  cout << "\nWriting geometry in specified format file...\n";
  mesh_atom.BuildIOFile(program,version,title,input,
			IOBaseT::kEnSight);
  

}


