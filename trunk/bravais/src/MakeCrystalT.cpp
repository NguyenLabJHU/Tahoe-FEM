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
  cout << "\n\nDimension: " << nsd << "\n";

  int whichunit;
  iArrayT cel(nsd);
  dArray2DT len(nsd,2);

  in >> whichunit;
  if(whichunit==1)
    {
      if(nsd==2) 
	{
	  in >> len(0,0) >> len(0,1);
	  in >> len(1,0) >> len(1,1);
	  cout << "Length read: [" << len(0,0) << "," << len(0,1) << "]\n";
	  cout << "             [" << len(1,0) << "," << len(1,1) << "]\n";
	}
      else if (nsd==3)
	{
	  in >> len(0,0) >>len(0,1);
	  in >> len(1,0) >>len(1,1);
	  in >> len(2,0) >>len(2,1);
	  cout << "Length read: [" << len(0,0) << "," << len(0,1) << "]\n";
	  cout << "             [" << len(1,0) << "," << len(1,1) << "]\n";
	  cout << "             [" << len(2,0) << "," << len(2,1) << "]\n";
	}
    }
  else
    {
      if(nsd==2) 
	{
	  in >> cel[0] >> cel[1];
	  cout << "Number of cells read: " << cel[0] << "  " << cel[1] <<  "\n";
	}
      else if (nsd==3)
	{
	  in >> cel[0] >> cel[1] >> cel[2];      
	  cout << "Number of cells read: " << cel[0] << "  " << cel[1] << "  " <<  cel[2] <<  "\n";
	}
    }
  StringT latticetype;
  in >> latticetype;
  cout << "Lattice Type: " << latticetype <<  "\n";
  
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

  StringT input = "example";

  //output format
  int intformat = 0;
  in >> intformat;
  IOBaseT::FileTypeT kformat = IOBaseT::int_to_FileTypeT(intformat);
  cout << "Output Format= " << kformat << "\n";

  iArrayT per(nsd);
  if (kformat == 12)
    {
      if(nsd == 2)
	{
	  in >> per[0] >> per[1] ;
	  cout << "Periodic conditions:\n";
	  cout << per[0] << "  " << per[1] << "\n";
	}
      else if(nsd == 3)
	{
	  in >> per[0] >> per[1] >> per[2];
	  cout << "Periodic conditions:\n";
	  cout << per[0] << "  " 
	       << per[1] << "  " 
	       << per[2] << "\n";
	}
    }


  //rotation
  double angle = 0.0;
  dArray2DT mat_rot(nsd,nsd);
  mat_rot = 0.0;

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
	  in >> mat_rot(0,0) >> mat_rot(1,0) >> mat_rot(2,0);
	  in >> mat_rot(0,1) >> mat_rot(1,1) >> mat_rot(2,1);
	  in >> mat_rot(0,2) >> mat_rot(1,2) >> mat_rot(2,2); 

	  cout << "Rotation Matrix:\n";
	  cout << mat_rot(0,0) << "  " <<  mat_rot(1,0) << "  " << mat_rot(2,0) << "\n";
	  cout << mat_rot(0,1) << "  " <<  mat_rot(1,1) << "  " << mat_rot(2,1) << "\n";
	  cout << mat_rot(0,2) << "  " <<  mat_rot(1,2) << "  " << mat_rot(2,2) << "\n";
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
		     shape,whichunit,len,cel,mat_rot,angle);

  StringT program = "bravais";
  StringT version = "v1.0";
  StringT title = "Lattice for Atoms";

  int nb_atoms;
  cout << "\nCreating mesh of atom...\n";
  nb_atoms = mesh_atom.CreateMeshAtom();
  cout << nb_atoms << " atoms in mesh\n";

  cout << "Final Lengths\n";
  if(nsd==2) 
    {
      cout << "[" << mesh_atom.Length()(0,0) << "," << mesh_atom.Length()(0,1) << "]\n";
      cout << "[" << mesh_atom.Length()(1,0) << "," << mesh_atom.Length()(1,1) << "]\n";
    }
  else if (nsd==3)
    {
      cout << "[" << mesh_atom.Length()(0,0) << "," << mesh_atom.Length()(0,1) << "]\n";
      cout << "[" << mesh_atom.Length()(1,0) << "," << mesh_atom.Length()(1,1) << "]\n";
      cout << "[" << mesh_atom.Length()(2,0) << "," << mesh_atom.Length()(2,1) << "]\n";
    }
  

  cout << "\nThe total volume of the domain is " 
       << mesh_atom.Volume_of_Mesh()
       << " in Angstroms^3 \n";

  dArray2DT coords;
  coords.Dimension(nb_atoms,nsd);

  coords = *(mesh_atom.ReturnCoordinates());

  if(nb_atoms <= 50) 
    {
      cout << "\nCoordinates:\n";
      if(nsd==2)
	{
	  for (int j=0; j<nb_atoms; j++) 
	    cout << j << " " << coords(j)[0] <<  "  " << coords(j)[1] <<  "\n";
	}
      else 
	{
	  for (int j=0; j<nb_atoms; j++) 
	    cout << j << " " <<  coords(j)[0] <<  "  " << coords(j)[1] <<  "  " << coords(j)[2] << "\n";
	}
    }

  cout << "\nWriting geometry in specified format file...\n";
  mesh_atom.BuildIOFile(program,version,title,input,kformat,per);
  

}


