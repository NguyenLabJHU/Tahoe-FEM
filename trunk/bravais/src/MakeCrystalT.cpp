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

  // Get name of input "data" file
  StringT inputfile;
  cout << "Name of input file?" << "\n";
  cin >> inputfile;
 
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

  double angle = 0.0;
  dArrayT rot_vec(nsd);
  rot_vec = 0.0;

  StringT input = "example";

  StringT misc;
  in >> misc;
 
  while(misc!="#")
  { 
     if(misc=="ROT") 
     {
        if(nsd==2) 
        {
           in >> angle;
           cout << "Rotation Angle: " << angle << "\n";
        }
        else if (nsd==3)
        {
	  in >> rot_vec[0] >> rot_vec[1] >> rot_vec[2] >> angle;
           cout << "Rotation Vector: " << rot_vec[0] 
		<< "  " <<  rot_vec[1] << "  " << rot_vec[2] 
		<< "\n";
	  cout  << "Angle: " << angle  << "\n";
        }
     }
     else if(misc=="OUTPUT") 
     {
        in >> input;
     }

     in >> misc;
  }

  in.close();
  
  //Define Mesh

  dArrayT latticeparameter;
  latticeparameter.Dimension(nsd);

  for (int i=0; i<nsd; i++)   
    latticeparameter[i] = alat;


  MeshAtom mesh_atom(latticetype,nsd,b,latticeparameter,
		     shape,whichunit,len_cell,rot_vec,angle);

  StringT program = "bravais";
  StringT version = "v1.0";
  StringT title = "Lattice for Atoms";

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
  /*
  if(nsd==2)
    {
      for (int j=0; j<nb_atoms; j++) 
	cout << coords(j,0) <<  "  " 
	     << coords(j,1) <<  "\n";
    }
  else 
    {
      for (int j=0; j<nb_atoms; j++) 
	cout << coords(j,0) <<  "  " 
	     << coords(j,1) <<  "  " 
	     << coords(j,2) << "\n";
    }
  */

  cout << "\nWriting geometry in specified format file...\n";
  mesh_atom.BuildIOFile(program,version,title,input,
			IOBaseT::kEnSight);
  

}


