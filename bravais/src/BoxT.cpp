/* $Id: BoxT.cpp,v 1.6 2002-09-09 23:10:29 saubry Exp $ */
#include "BoxT.h"
#include "VolumeT.h"

#include <iostream>
#include <fstream>

#include "ExceptionCodes.h"
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
  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();
  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vB = pcl->GetBasis();
  const dArray2DT& vA = pcl->GetAxis();

  int natoms,a = 0;
  int temp_nat=0;
  dArray2DT temp_atom;
  double x,y,z;
  double eps = 1.e-6;

  // Define a slightly shorter box
  double l00,l01,l10,l11;
  if(length(0,0) >= 0.0) l00 = length(0,0) + eps; else l00 = length(0,0) - eps;
  if(length(0,1) >= 0.0) l01 = length(0,1) + eps; else l01 = length(0,1) - eps;
  if(length(1,0) >= 0.0) l10 = length(1,0) + eps; else l10 = length(1,0) - eps;
  if(length(1,1) >= 0.0) l11 = length(1,1) + eps; else l11 = length(1,1) - eps;

  if (nlsd==2) temp_nat = 16*nuca*ncells[0]*ncells[1];
  if (nlsd==3) temp_nat = 64*nuca*ncells[0]*ncells[1]*ncells[2];

  temp_atom.Dimension(temp_nat,nlsd);

  if (nlsd==2) 
    {
      // Rotate coordinates
      for (int p=-2*ncells[1];p<2*ncells[1];p++) 
      	for (int q=-2*ncells[0];q<2*ncells[0];q++)    
	  {
	    dArrayT c(nlsd);
	    c[0] = (double)q; c[1] = (double)p;
	    for (int m=0;m<nuca;m++) 
	      {
		if (a >= temp_nat) {throw eSizeMismatch;}
		
		x = length(0,0);
		y = length(1,0);
		
		for (int k=0;k<nlsd;k++) 
		  {
		    x += (c[k] + vB(k,m))*vA(k,0);
		    y += (c[k] + vB(k,m))*vA(k,1);
		  }

		//		cout << x << "  " << y << "\n";

		if(x >= l00 && x <= l01 && y >= l10 && y <= l11 )
		  {
		    temp_atom(a)[0] = x;
		    temp_atom(a)[1] = y;
		    a++;
		  }
	      }
	  }
    }
  else if (nlsd==3) 
    {
      double l20,l21;
      if(length(2,0) >= 0.0) l20 = length(2,0) + eps; else l20 = length(2,0) - eps;
      if(length(2,1) >= 0.0) l21 = length(2,1) + eps; else l21 = length(2,1) - eps;

      // Rotate coordinates
      for (int p=-2*ncells[2];p<2*ncells[2];p++) 
	for (int q=-2*ncells[1];q<2*ncells[1];q++)    
	  for (int r=-2*ncells[0];r<2*ncells[0];r++) 
	    {
	      dArrayT c(nlsd);
	      c[0] = (double)r; c[1] = (double)q; c[2] = (double)p;
	      for (int m=0;m<nuca;m++) 
		{
		  if (a > temp_nat) {throw eSizeMismatch;}

		  x = length(0,0);
		  y = length(1,0);
		  z = length(2,0);
		  
		  for (int k=0;k<nlsd;k++) 
		    {
		      x += (c[k] + vB(k,m))*vA(k,0);
		      y += (c[k] + vB(k,m))*vA(k,1);
		      z += (c[k] + vB(k,m))*vA(k,2);
		    }
		  
		  if(x >= l00 && x <= l01 && y >= l10 && y <= l11 &&
		     z >= l20 && z <= l21)
		    {
		      temp_atom(a)[0] = x;
		      temp_atom(a)[1] = y;
		      temp_atom(a)[2] = z;
		      a++;                     
		    }
		}
	    }
    }

  // Get atoms coordinates
  nATOMS = a;
  atom_ID.Dimension(nATOMS);
  atom_coord.Dimension(nATOMS,nlsd);
  atom_connectivities.Dimension(nATOMS,1);

  for(int m=0; m < nATOMS ; m++) 
    {
      for (int k=0;k< nlsd;k++)
	atom_coord(m)[k] = temp_atom(m)[k];
    }
  
  atom_names = "Box";
  for (int p=0;p<nATOMS;p++)
    {
      atom_ID[p] = p;
      atom_connectivities(p)[0] = p;
    }

  // Update lengths
  length = ComputeMinMax();

  //Calculate volume here
  switch(nlsd) {
  case 2:
    volume = (length(0,1)-length(0,0))*(length(1,1)-length(1,0));
    break;
  case 3:
    volume = (length(0,1)-length(0,0))*         
             (length(1,1)-length(1,0))*
             (length(2,1)-length(2,0));
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
    {
      if(minmax(i,0) > minmax(i,1)) 
	{
	  double temp = minmax(i,0);
	  minmax(i,0) = minmax(i,1);
	  minmax(i,1) = temp;
	}
    }

  return minmax;  
}
