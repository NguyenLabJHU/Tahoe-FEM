/* $Id: BoxT.cpp,v 1.4 2002-07-25 23:48:07 saubry Exp $ */
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

BoxT::BoxT(int dim, int whichunit, dArrayT len_cel, 
	   dArrayT lattice_parameter) : VolumeT(dim) 
{
  nSD = dim;
  length.Dimension(nSD);
  ncells.Dimension(nSD);

  for(int i=0;i<nSD;i++)
    {
      length[i] = len_cel[i]*lattice_parameter[i];
      ncells[i] = static_cast<int>(len_cel[i]);
    }
}


BoxT::BoxT(const BoxT& source) : VolumeT(source.nSD) 
{
  ncells.Dimension(source.nSD);
  ncells = source.ncells;

  length.Dimension(source.nSD);
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

void BoxT::CalculateVolume() 
{
  switch(nSD) {
  case 2:
    volume = length[0]*length[1];
    break;
  case 3:
    volume = length[0]*length[1]*length[2];
    break;
  }
}

void BoxT::CreateLattice(CrystalLatticeT* pcl) 
{
  int nuca = pcl->GetNUCA();
  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vB = pcl->GetBasis();
  int p,q,r,a = 0;
  
  
  cout << "nLSD=" << nSD << "\n";
  cout << "nUCA=" << nuca << "\n";
  
  if (nSD==2) nATOMS = nuca*ncells[0]*ncells[1];
  if (nSD==3) nATOMS = nuca*ncells[0]*ncells[1]*ncells[2];
  
  atom_ID.Dimension(nATOMS);
  atom_coord.Dimension(nATOMS,nSD);
  atom_connectivities.Dimension(nATOMS,1);
  
  if (nSD==2) 
    {
      for (p=-(ncells[1]/2);p<((ncells[1]+1)/2);p++) 
	{
	  for (q=-(ncells[0]/2);q<((ncells[0]+1)/2);q++) 
	    {
	      for (int m=0;m<nuca;m++) 
		{
		  
		  if (a > nATOMS) {cout << "nATOMS wrong";throw eSizeMismatch;}
		  
		  atom_coord(a)[0] = ((double) q + vB(0,m))*vLP[0];
		  atom_coord(a)[1] = ((double) p + vB(1,m))*vLP[1];
		  
		  a++;
    }
	    }
	}
    }
  else if (nSD==3) 
    {
      for (p=-((ncells[2])/2);p<((ncells[2]+1)/2);p++) 
	{
	  for (q=-((ncells[1])/2);q<((ncells[1]+1)/2);q++) 
	    {
	      for (r=-((ncells[0])/2);r<((ncells[0]+1)/2);r++) 
		{
		  for (int m=0;m<nuca;m++) 
		    {
		      
		      if (a > nATOMS) {cout << "nATOMS wrong";throw eSizeMismatch;}
		      
		      atom_coord(a)[0] = ((double) r + vB(0,m))*vLP[0];
		      atom_coord(a)[1] = ((double) q + vB(1,m))*vLP[1];
		      atom_coord(a)[2] = ((double) p + vB(2,m))*vLP[2];
		      
		      a++;                     }
		  
		}
	    }
	}
    }

  // Rotate coordinates:
  atom_coord = pcl->BasisRotation(atom_coord);
  
  
  atom_names = "Box";
  for (p=0;p<nATOMS;p++)
    {
      atom_ID[p] = p;
      atom_connectivities(p)[0] = p;
    }
  
}
