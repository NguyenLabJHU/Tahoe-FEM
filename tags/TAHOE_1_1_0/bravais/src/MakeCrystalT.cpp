// DEVELOPMENT
/* Build a mesh of atoms with a format-independent output
 * (saubry, Tue Jul 23 15:31 UTC 2002)
*/


#include "MakeCrystalT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ExceptionCodes.h"

#include "VolumeT.h"
#include "BoxT.h"

#include "CrystalLatticeT.h"
#include "CUBT.h"
#include "FCCT.h"
#include "BCCT.h"
#include "DIAT.h"
#include "HEXT.h"

#include "OutputSetT.h"
#include "OutputBaseT.h"
#include "OutPutLatticeT.h"

#include "MeshAtom.h"

void MakeCrystalT::Run() 
{
  // Get name of input data file
  StringT inputfile;
  cin >> inputfile;

  ifstreamT in('%');
  in.open(inputfile);

  // Create output file
  StringT outputfile;
  outputfile = inputfile.Root();
  outputfile.Append(".out");
  ofstreamT out;
  out.open(outputfile);
 
  // Start to read data
  int nsd=0; 
  in >> nsd;
  out << "\n\nDimension: " << nsd << "\n";
 
  int whichunit;
  iArrayT cel(nsd);
  dArray2DT len(nsd,2);
 
  cel = 0;
  len = 0.0;

  in >> whichunit;
  if(whichunit==1)
    {
      if(nsd==2) 
	{
	  in >> len(0,0) >> len(0,1);
	  in >> len(1,0) >> len(1,1);
	  out << "Length read: [" << len(0,0) << "," << len(0,1) << "]\n";
	  out << "             [" << len(1,0) << "," << len(1,1) << "]\n";
	}
      else if (nsd==3)
	{
	  in >> len(0,0) >>len(0,1);
	  in >> len(1,0) >>len(1,1);
	  in >> len(2,0) >>len(2,1);
	  out << "Length read: [" << len(0,0) << "," << len(0,1) << "]\n";
	  out << "             [" << len(1,0) << "," << len(1,1) << "]\n";
	  out << "             [" << len(2,0) << "," << len(2,1) << "]\n";
	}
    }
  else
    {
      if(nsd==2) 
	{
	  in >> cel[0] >> cel[1];
	  out << "Number of cells read: " << cel[0] << "  " << cel[1] <<  "\n";
	}
      else if (nsd==3)
	{
	  in >> cel[0] >> cel[1] >> cel[2];      
	  out << "Number of cells read: " << cel[0] << "  " << cel[1] << "  " <<  cel[2] <<  "\n";
	}
    }

  StringT latticetype;
  in >> latticetype;
  out << "Lattice Type: " << latticetype <<  "\n";
  int b=0;
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
  else if (latticetype=="HEX") 
    {
      if(nsd == 2) b=2;
      if(nsd == 3) b=3;
    }
  else if (latticetype=="CUB") 
    {
      if(nsd == 2) b=1;
      if(nsd == 3) b=1;
    }
  else 
    {
      out << "Lattice type has to be CUB, FCC, BCC, DIA or HEX...\n";
      throw eBadInputValue;
    }
	
  dArrayT alat(nsd);
  if(nsd == 2)
    {
      in >> alat[0] >> alat[1];
      out << "Lattice parameter: " 
	   << alat[0] << "  " 
	   << alat[1] << "\n";
    }
  else
    {
      in >> alat[0] >> alat[1] >> alat[2];
      out << "Lattice parameter: " 
	   << alat[0] << "  " 
	   << alat[1] << "  " 
	   << alat[2] << "\n";
    }

  StringT shape;
  in >> shape;
  out << "Shape of the domain:" << shape <<  "\n";

  //read output format
  int intformat = 0;
  in >> intformat;

  IOBaseT::FileTypeT kformat = IOBaseT::int_to_FileTypeT(intformat);
  out << "Output Format: " << kformat << "\n";
 
  //Set Defaults on periodicity, rotation and output filenames

  //periodicity -- default value is 0 (non-periodic boundary conditions)
  iArrayT per(nsd);
  per = 0;

  //rotation -- default value is zero rotation
  double angle = 0.0;
  dArray2DT mat_rot(nsd,nsd);
  mat_rot = 0.0;
  
  for(int i=0;i<nsd;i++) 
    mat_rot(i,i) = 1.0;

  //output filename prefix
  StringT input = "example";

  //Override Defaults using miscellaneous input arguments
  //or use the "#" symbol to end the reading of the data file

  StringT misc;
  in >> misc;

  int irot= 0;

  iArrayT isort(nsd);
  isort = 0;

  while (misc!="#")
  {
    if (misc=="PERIODICITY")
    {
      if(nsd == 2)
      {
        in >> per[0] >> per[1] ;
        out << "Periodic conditions:\n";
        out << per[0] << "  " << per[1] << "\n";
      }
      else if(nsd == 3)
      {
        in >> per[0] >> per[1] >> per[2];
        out << "Periodic conditions:\n";
        out << per[0] << "  "  << per[1] << "  " << per[2] << "\n";
      }
    }
    else if (misc=="ROTATION")
    {
      in >> irot;
      if (irot == 0) out << "Rotation atoms in box\n";
      if (irot == 1) out << "Rotation box of atoms\n";

      if(nsd==2) 
      {
        in >> angle;
        out << "Rotation Angle: " << angle << "\n";
      }
      else if (nsd==3)
      {
	// Read in the three vectors written in line, 
	// set them up in a matrix as column vectors.
        in >> mat_rot(0,0) >> mat_rot(1,0) >> mat_rot(2,0);
        in >> mat_rot(0,1) >> mat_rot(1,1) >> mat_rot(2,1);
        in >> mat_rot(0,2) >> mat_rot(1,2) >> mat_rot(2,2); 
      
        out << "Rotation Matrix:\n";
        out << mat_rot(0,0) << "  " <<  mat_rot(0,1) << "  " << mat_rot(0,2) << "\n";
        out << mat_rot(1,0) << "  " <<  mat_rot(1,1) << "  " << mat_rot(1,2) << "\n";
        out << mat_rot(2,0) << "  " <<  mat_rot(2,1) << "  " << mat_rot(2,2) << "\n";
      }
    }
    else if (misc=="SORT")
      {
	if(nsd == 2) 
	  {
	    in >> isort[0] >> isort[1];
	    out << "Sorting criteria: ";
	    out << "[" << isort[0] << " " << isort[1] << "]\n";
	  }
	if(nsd == 3) 
	  {
	    in >> isort[0] >> isort[1] >> isort[2];
	    out << "Sorting criteria: ";
	    out << "[" << isort[0] << " " << isort[1] 
		 << " " << isort[2] << "]\n";
	  }
      }
    else if (misc=="OUTPUT")
    {
      in >> input;
      out << "Output file root: " << input << "\n";
    } 

    in >> misc;
  }

  // Close data file
  in.close();

  //Define Mesh

  MeshAtom mesh_atom(latticetype,nsd,b,alat,
		     shape,whichunit,len,cel,irot,mat_rot,
		     angle,isort);

  StringT program = "bravais";
  StringT version = "v1.0";
  StringT title = "Lattice for Atoms";

  int nb_atoms;
  out << "\nCreating mesh of atom...\n";
  nb_atoms = mesh_atom.CreateMeshAtom();
  out << nb_atoms << " atoms in mesh\n";

  out << "\nActual Length of Atoms\n";
  if(nsd==2) 
    {
      out << "[" << mesh_atom.Length()(0,0) << "," << mesh_atom.Length()(0,1) << "]\n";
      out << "[" << mesh_atom.Length()(1,0) << "," << mesh_atom.Length()(1,1) << "]\n";
    }
  else if (nsd==3)
    {
      out << "[" << mesh_atom.Length()(0,0) << "," << mesh_atom.Length()(0,1) << "]\n";
      out << "[" << mesh_atom.Length()(1,0) << "," << mesh_atom.Length()(1,1) << "]\n";
      out << "[" << mesh_atom.Length()(2,0) << "," << mesh_atom.Length()(2,1) << "]\n";
    }

  out << "\nThe total volume of the domain is " 
       << mesh_atom.Volume_of_Mesh()
       << " in Angstroms^3 \n";

  dArray2DT coords;
  coords.Dimension(nb_atoms,nsd);

  coords = *(mesh_atom.ReturnCoordinates());

  if(nb_atoms <= 50) 
    {
      out << "\nCoordinates:\n";
      if(nsd==2)
	{
	  for (int j=0; j<nb_atoms; j++) 
	    out << coords(j)[0] <<  "  " << coords(j)[1]<<  "\n";
	}
      else 
	{
	  for (int j=0; j<nb_atoms; j++) 
	    out << coords(j)[0] <<  "  " << coords(j)[1] <<  "  " << coords(j)[2] << "\n";
	}
      
    }


  //Output in a file
  /*
  StringT outputcoord = "coord.txt";
  ofstreamT out_coord;
  out_coord.open(outputcoord);
  
  if(nsd==2)
    {
      for (int j=0; j<nb_atoms; j++) 
	out_coord << coords(j)[0] << "  " << coords(j)[1] <<  "\n";
    }
  else 
    {
      for (int j=0; j<nb_atoms; j++) 
	out_coord << coords(j)[0] << "  " << coords(j)[1] 
		  << "  " <<coords(j)[2] <<  "\n";
    }

  out_coord.close();

  */

  out << "\nWriting geometry in specified format file...\n";
  mesh_atom.BuildIOFile(program,version,title,input,kformat,per);

}
