/* $Id: BoxT.cpp,v 1.5 2002-08-02 02:07:49 saubry Exp $ */
#include "BoxT.h"
#include "VolumeT.h"

#include <iostream>
#include <fstream>

#include "ifstreamT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

BoxT::BoxT(int dim, dArray2DT len,
	   dArrayT lattice_parameter) : VolumeT(dim) 
{
  nSD = dim;
  length.Dimension(nSD,2);
  ncells.Dimension(nSD);

  for(int i=0;i<nSD;i++)
    {
      double dist = len(i,1)-len(i,0);
      if(dist < lattice_parameter[i])
	{
	  cout << "Lengths too small to define atoms lattice\n";
	  cout << "Lengths have to be greater than lattice parameter\n";
	  cout << "i=" << i << "dist=" << dist;
	  throw eBadInputValue;
	}
      
      ncells[i] = static_cast<int>(dist/lattice_parameter[i]);
    }

  length = len;
}

BoxT::BoxT(int dim, iArrayT cel,
	   dArrayT lattice_parameter) : VolumeT(dim) 
{
  nSD = dim;
  length.Dimension(nSD,2);
  ncells.Dimension(nSD);

  for(int i=0;i<nSD;i++)
      ncells[i] = cel[i];

  for(int i=0;i<nSD;i++)
    {
      length(i,0) = -cel[i]*lattice_parameter[i]*0.5;
      length(i,1) =  cel[i]*lattice_parameter[i]*0.5;
    }
}


BoxT::BoxT(const BoxT& source) : VolumeT(source.nSD) 
{
  ncells.Dimension(source.nSD);
  ncells = source.ncells;

  length.Dimension(source.nSD,2);
  length = source.length;

  volume = source.volume;

  atom_names = source.atom_names;

  atom_ID.Dimension(source.nATOMS);
  atom_ID = source.atom_ID;

  atom_coord.Dimension(source.nATOMS,source.nSD);
  atom_coord = source.atom_coord;


  atom_connectivities.Dimension(source.nATOMS,source.nSD);
  atom_connectivities = source.atom_connectivities;
}

void BoxT::CreateLattice(CrystalLatticeT* pcl) 
{
  int nuca = pcl->GetNUCA();
  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vB = pcl->GetBasis();
  int p,q,r,a = 0;
  
  if (nSD==2) nATOMS = nuca*ncells[0]*ncells[1];
  if (nSD==3) nATOMS = nuca*ncells[0]*ncells[1]*ncells[2];

  if(nATOMS == 0) 
    {
      cout << "number of atoms is zero\n";
      throw eBadInputValue;
    }
  
  atom_ID.Dimension(nATOMS);
  atom_coord.Dimension(nATOMS,nSD);
  atom_connectivities.Dimension(nATOMS,1);

  if (nSD==2) 
    {
      for (p=0;p<ncells[1];p++) 
	for (q=0;q<ncells[0];q++)    
	  {
	    for (int m=0;m<nuca;m++) 
	      {
		if (a > nATOMS) {cout << "nATOMS wrong\n";throw eSizeMismatch;}
		
		atom_coord(a)[0] = length(0,0) + ((double) q + vB(0,m))*vLP[0];
		atom_coord(a)[1] = length(1,0) + ((double) p + vB(1,m))*vLP[1];
		
		a++;
	      }
	  }
    }
  else if (nSD==3) 
    {
      for (p=0;p<ncells[2];p++) 
	for (q=0;q<ncells[1];q++)    
	  for (r=0;r<ncells[0];r++) 		
	    for (int m=0;m<nuca;m++) 
	      {
		if (a > nATOMS) {cout << "nATOMS wrong\n";throw eSizeMismatch;}
		
		atom_coord(a)[0] = length(0,0) + ((double) r + vB(0,m))*vLP[0];
		atom_coord(a)[1] = length(1,0) + ((double) q + vB(1,m))*vLP[1];
		atom_coord(a)[2] = length(2,0) + ((double) p + vB(2,m))*vLP[2];
		
		a++;                     
	      }
    }
  
  atom_names = "Box";
  for (p=0;p<nATOMS;p++)
    {
      atom_ID[p] = p;
      atom_connectivities(p)[0] = p;
    }


  // Rotate coordinates:
  atom_coord = pcl->BasisRotation(atom_coord);

  // Rotate Axis
  cout << "Previous Lengths\n";
  if(nSD==2) 
    {
      cout << "[" << length(0,0) << "," << length(0,1) << "]\n";
      cout << "[" << length(1,0) << "," << length(1,1) << "]\n";
    }
  else if (nSD==3)
    {
      cout << "[" << length(0,0) << "," << length(0,1) << "]\n";
      cout << "[" << length(1,0) << "," << length(1,1) << "]\n";
      cout << "[" << length(2,0) << "," << length(2,1) << "]\n";
    }

  length = ComputeMinMax();

  //Calculate volume here
  switch(nSD) {
  case 2:
    volume = (length(0,1)-length(0,0))*(length(1,1)-length(1,0));
    break;
  case 3:
    volume = (length(0,1)-length(0,0))*(length(1,1)-length(1,0))*(length(2,1)-length(2,0));
    break;
  }
}

void BoxT::CalculateBounds(iArrayT per,CrystalLatticeT* pcl)
{
  const dArrayT& vLP = pcl->GetLatticeParameters();
  //dArray2DT MinMax(nSD,2);

  //MinMax = ComputeMinMax();

  atom_bounds.Dimension(nSD,2);

  for (int i=0; i < nSD; i++)
    {
      if (per[i]==0) 
	{
	  // non-periodic conditions
	  atom_bounds(i,0) = -10000.;
	  atom_bounds(i,1) =  10000.;
	}
      else if (per[i]==1)
	{
	  // periodic conditions
	  //	  atom_bounds(i,0) = MinMax(i)[0];
	  //	  atom_bounds(i,1) = MinMax(i)[1] + 0.5*vLP[1];

	  	  atom_bounds(i,0) = length(i)[0];
	  	  atom_bounds(i,1) = length(i)[1] + 0.5*vLP[1];
	}
      else
	throw eBadInputValue;
    }  
}


void BoxT::CalculateType()
{
  atom_types.Dimension(nATOMS);

  for (int i=0; i < nATOMS; i++)
    atom_types[i] = 1; 
}



//////////////////// PRIVATE //////////////////////////////////

dArray2DT BoxT::ComputeMinMax()
{
  dArray2DT minmax(nSD,2);

  for (int i=0; i < nSD; i++)
    {
      minmax(i,0) = atom_coord(0)[i]; // min
      minmax(i,1) = atom_coord(0)[i]; //max
    }
  
  for (int i=0; i < nSD; i++)
      for (int j=1; j < nATOMS; j++)
	{
	  if(atom_coord(j)[i] < minmax(i,0)) minmax(i,0)=atom_coord(j)[i];
	  if(atom_coord(j)[i] > minmax(i,1)) minmax(i,1)=atom_coord(j)[i];	  
	}

  for (int i=0; i < nSD; i++)
    if(minmax(i,0) > minmax(i,1)) 
      {
	double temp = minmax(i,0);
	minmax(i,0) = minmax(i,1);
	minmax(i,1) = temp;
      }

  return minmax;  
}
