// DEVELOPMENT
/* $Id: BoxT.cpp,v 1.34 2003-08-01 22:54:39 saubry Exp $ */
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
	   dArrayT lattice_parameter,
	   iArrayT which_sort, StringT slt,
           iArrayT per) : VolumeT(dim) 
{
  nSD = dim;
  length.Dimension(nSD,2);
  ncells.Dimension(nSD);

  WhichSort.Dimension(nSD);
  WhichSort = which_sort;

  pbc.Dimension(nSD);
  pbc = per;

  sLATTYPE = slt;

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

  if (sLATTYPE == "CORUN")
    {
      double dist[3];
      dist[0] = 3.0*ncells[0]*lattice_parameter[0];
      dist[1] = sqrt(3.0)*ncells[1]*lattice_parameter[1];
      dist[2] = ncells[2]*lattice_parameter[2];
      for(int i=0;i<nSD;i++)
        {
          dist[i] = 0.5*dist[i];
          length(i,0) = -dist[i];
          length(i,1) = length(i,0) + (dist[i] - length(i,0));
        }
    }
  else if (sLATTYPE == "HEX")
    {
      double dist[3]; 
      dist[0] = ncells[0]*lattice_parameter[0];
      dist[1] = sqrt(3.0)*ncells[1]*lattice_parameter[1];
      if (nSD==3) dist[2] = ncells[2]*lattice_parameter[2];
      for(int i=0;i<nSD;i++)
        {
          dist[i] = 0.5*dist[i];                           
          length(i,0) = -dist[i];
          length(i,1) = length(i,0) + (dist[i] - length(i,0)); 
        }
    }    
  else
    {
      for(int i=0;i<nSD;i++)
        {
          double dist = ncells[i]*lattice_parameter[i]*0.5;
          length(i,0) = -dist;
          length(i,1) = length(i,0) + (dist - length(i,0));
        }
     }
}

BoxT::BoxT(int dim, iArrayT cel,
	   dArrayT lattice_parameter,
	   iArrayT which_sort, StringT slt,
           iArrayT per) : VolumeT(dim) 
{
  nSD = dim;
  length.Dimension(nSD,2);
  ncells.Dimension(nSD);

  WhichSort.Dimension(nSD);
  WhichSort = which_sort;

  pbc.Dimension(nSD);
  pbc = per;

  sLATTYPE = slt;

  for(int i=0;i<nSD;i++)
      ncells[i] = cel[i];

  if (sLATTYPE == "CORUN")
    {
      double dist[3];
      dist[0] = 3.0*ncells[0]*lattice_parameter[0];
      dist[1] = sqrt(3.0)*ncells[1]*lattice_parameter[1];
      dist[2] = ncells[2]*lattice_parameter[2];
      for(int i=0;i<nSD;i++)
        {
          dist[i] = 0.5*dist[i];
          length(i,0) = -dist[i];
          length(i,1) = length(i,0) + (dist[i] - length(i,0));
        }
    }
  else if (sLATTYPE == "HEX")
    {
      double dist[3];
      dist[0] = ncells[0]*lattice_parameter[0];
      dist[1] = sqrt(3.0)*ncells[1]*lattice_parameter[1];
      if (nSD==3) dist[2] = ncells[2]*lattice_parameter[2];
      for(int i=0;i<nSD;i++)
        {
          dist[i] = 0.5*dist[i];
          length(i,0) = -dist[i];
          length(i,1) = length(i,0) + (dist[i] - length(i,0));
        }
    }
  else
    {
      for(int i=0;i<nSD;i++)
        {
          double dist = ncells[i]*lattice_parameter[i]*0.5;
          length(i,0) = -dist;
          length(i,1) = length(i,0) + (dist - length(i,0));
        }
     }
}

BoxT::BoxT(const BoxT& source) : VolumeT(source.nSD) 
{
  ncells.Dimension(source.nSD);
  ncells = source.ncells;

  length.Dimension(source.nSD,2);
  length = source.length;

  volume = source.volume;

  WhichSort.Dimension(nSD);
  WhichSort = source.WhichSort;

  pbc.Dimension(source.nSD);
  pbc = source.pbc;

  atom_names = source.atom_names;

  atom_ID.Dimension(source.nATOMS);
  atom_ID = source.atom_ID;

  atom_coord.Dimension(source.nATOMS,source.nSD);
  atom_coord = source.atom_coord;

  atom_types.Dimension(source.nATOMS);
  atom_types = source.atom_types;

  atom_connectivities.Dimension(source.nATOMS,source.nSD);
  atom_connectivities = source.atom_connectivities;
}


void BoxT::CreateLattice(CrystalLatticeT* pcl) 
{
  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();
  int ntype = pcl->GetNTYPE();

  int natoms=0;
  int temp_nat=0;
  dArray2DT temp_atom;
  iArrayT temp_type;

  if(pcl->GetRotMeth() == 0) 
    {
      if (nlsd==2) temp_nat = 16*nuca*ncells[0]*ncells[1];
      if (nlsd==3) temp_nat = 64*nuca*ncells[0]*ncells[1]*ncells[2];
    }
  else
    {
      if (nlsd==2) temp_nat =  8*nuca*ncells[0]*ncells[1];
      if (nlsd==3) temp_nat = 16*nuca*ncells[0]*ncells[1]*ncells[2];
    }
  temp_atom.Dimension(temp_nat,nlsd);
  temp_type.Dimension(temp_nat);

  if(pcl->GetRotMeth() == 0)
    nATOMS = RotateAtomInBox(pcl,&temp_atom,&temp_type,temp_nat);
  else
    nATOMS = RotateBoxOfAtom(pcl,&temp_atom,&temp_type,temp_nat);


  // Get atoms coordinates
  atom_ID.Dimension(nATOMS);
  atom_coord.Dimension(nATOMS,nlsd);
  atom_types.Dimension(nATOMS);
  atom_connectivities.Dimension(nATOMS,1);

  for(int m=0; m < nATOMS ; m++) 
    for (int k=0;k< nlsd;k++)
      {
	atom_coord(m)[k] = temp_atom(m)[k];
	atom_types[m] = temp_type[m];
      }
  
  atom_names = "Box";

  // Create connectivities and IDs
  for (int p=0;p<nATOMS;p++)
    {
      atom_ID[p] = p;
      atom_connectivities(p)[0] = p;
    }

  if(ntype > 2) 
    cout << "WARNING: ** nTypes == 2  maximum for ensight output ** \n"

  atom_array_ID.Dimension(ntype);
  atom_array_connect.Dimension(ntype);


  type1.Dimension(nATOMS,1);type1 = 0;
  type2.Dimension(nATOMS,1);type2 = 0;

  int n=0,m=0;
  for (int p=0;p<nATOMS;p++)
    {
      if ( atom_types[p] == 1) 
	{
	  type1(n)[0] = p;
	  n++;
	}
      else
	{
	  type2(m)[0] = p;
	  m++;
	}	
    }

  atom_array_ID[0] = "type1";
  if (ntype > 1) atom_array_ID[1] = "type2";

  if (n < nATOMS && n > 0) type1.Resize(n);
  if (m < nATOMS && m > 0) type2.Resize(m);
  atom_array_connect[0] = &type1;
  if (ntype > 1) atom_array_connect[1] = &type2;


  // Create parts
  atom_parts.Dimension(nATOMS);
  atom_parts = 1;

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

  // Sort Lattice 
  if(WhichSort  != 0) SortLattice(pcl);

}

void BoxT::SortLattice(CrystalLatticeT* pcl)
{
  int nlsd = pcl->GetNLSD();
  
  dArray2DT new_coord;
  new_coord.Dimension(atom_coord.MajorDim(),nlsd);   
  new_coord = 0.0;

  iArrayT new_type;
  new_type.Dimension(nATOMS);  
  new_type = 0;

  dArrayT x(atom_coord.MajorDim()); x = 0.0;
  dArrayT y(atom_coord.MajorDim()); y = 0.0;
  dArrayT z(atom_coord.MajorDim()); z = 0.0;

  iArrayT Map(atom_coord.MajorDim());
  Map = 0;

  for(int m=0; m < nATOMS ; m++) 
    {
      x[m] = atom_coord(m)[WhichSort[0]];
      y[m] = atom_coord(m)[WhichSort[1]];
      if (nlsd == 3) z[m] = atom_coord(m)[WhichSort[2]];
    }

  // Sort 1st criterium
  Map.SetValueToPosition();
  Map.SortAscending(x);

  for(int m=0; m < nATOMS ; m++) 
    {
      new_coord(m)[WhichSort[0]] = x[m];
      new_coord(m)[WhichSort[1]] = y[Map[m]];
      if (nlsd == 3)  new_coord(m)[WhichSort[2]] = z[Map[m]];
      new_type[m] = atom_types[Map[m]];
    } 

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;

  // Sort 2nd criterium
  new_coord = 0.0;
  new_type = 0;

  iArrayT Ind(atom_coord.MajorDim());
  Ind = 0;

  int is = 0;
  Ind[is] = 0;
  for(int m = 1; m < nATOMS ; m++) 
    {
      int ws = WhichSort[0];
      if( fabs(atom_coord(m)[ws] - atom_coord(m-1)[ws]) >  1.e-10 ) 
	{
	  is++;
	  Ind[is] = m ;
	}
    }  
  is++;
  Ind[is] = nATOMS;
  is++;

  x = 0;
  y = 0;
  z = 0;
  int p = 0;
  for(int n = 0; n < is-1; n++)
    {
      dArrayT aux(Ind[n+1]-Ind[n]);
      dArrayT aux2(Ind[n+1]-Ind[n]);
      aux = 0.0;
      aux2= 0.0;

      int isa = 0;
      for(int m = Ind[n]; m < Ind[n+1]; m++)
	{
	  aux[isa] = atom_coord(m)[WhichSort[1]];
	  if(nlsd == 3) aux2[isa]= atom_coord(m)[WhichSort[2]];
	  isa++;
	}

      iArrayT Map2(isa);
      Map2.SetValueToPosition();
      Map2.SortAscending(aux);

      for(int m = 0; m < isa; m++)
	{
	  y[p] = aux[m];
	  if(nlsd == 3) z[p] = aux2[Map2[m]];
	  new_type[p] = atom_types[Map2[m]];
	  p++;
	}
    }
  
  for(int m=0; m < nATOMS ; m++) 
    {
      new_coord(m)[WhichSort[0]] = atom_coord(m)[WhichSort[0]];
      new_coord(m)[WhichSort[1]] = y[m];
      if (nlsd == 3) new_coord(m)[WhichSort[2]] = z[m];
    } 

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;

  // Sort 3nd criterium
  if(nlsd == 3)
  {
    new_coord = 0.0;
    new_type = 0;

    iArrayT Ind(atom_coord.MajorDim());
    Ind = 0;

    int is = 0;
    Ind[is] = 0;
    for(int m = 1; m < nATOMS ; m++) 
      {
	int ws = WhichSort[1];
	if( fabs(atom_coord(m)[ws] - atom_coord(m-1)[ws]) >  1.e-10 ) 
	  {
	    is++;
	    Ind[is] = m ;
	  }
      }  
    is++;
    Ind[is] = nATOMS;
    is++;
    
    x = 0;
    y = 0;
    z = 0;
    int p = 0;
    for(int n = 0; n < is-1; n++)
      {
	dArrayT aux(Ind[n+1]-Ind[n]);
	aux = 0.0;

	int isa = 0;
	for(int m = Ind[n]; m < Ind[n+1]; m++)
	  {
	    aux[isa] = atom_coord(m)[WhichSort[2]];
	    isa++;
	  }
	
	iArrayT Map2(isa);
	Map2.SetValueToPosition();
	Map2.SortAscending(aux);
	
	for(int m = 0; m < isa; m++)
	  {
	    z[p] = aux[m];
	    new_type[p] = atom_types[m];
	    p++;
	  }
      }
    
    for(int m=0; m < nATOMS ; m++) 
      {
	new_coord(m)[WhichSort[0]] = atom_coord(m)[WhichSort[0]];
	new_coord(m)[WhichSort[1]] = atom_coord(m)[WhichSort[1]];
	new_coord(m)[WhichSort[2]] = z[m];
      } 
  }

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;

}

void BoxT::CalculateBounds(CrystalLatticeT* pcl)
{
  // const dArrayT& vLP = pcl->GetLatticeParameters();
  // const dArray2DT& vAX = pcl->GetAxis();

  atom_bounds.Dimension(nSD,2);

  for (int i=0; i < nSD; i++)
    {
      if (pbc[i]==0) 
	{
	  // non-periodic conditions
	  atom_bounds(i,0) = -10000.;
	  atom_bounds(i,1) =  10000.;
	}
      else if (pbc[i]==1)
	{
	  // periodic conditions
	  // atom_bounds(i,0) = length(i)[0];
	  // atom_bounds(i,1) = length(i)[1] + 0.5*vLP[1];
	  // atom_bounds(i,1) = length(i)[1] + 0.5*vAX(i,i);
          // cout << vAX(i,i) << "\n";
          atom_bounds(i,0) = length(i,0);
          atom_bounds(i,1) = length(i,1);
	}
      else
	throw eBadInputValue;
    }  
}


//////////////////// PRIVATE //////////////////////////////////

int BoxT::RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat)
{
  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();

  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vB = pcl->GetBasis();
  const dArray2DT& vA = pcl->GetAxis();
  const iArrayT& vT = pcl->GetType();

  double x,y,z;
  double eps = 1.e-6;

  int natom= 0;
  int type= 0;

  // Define a slightly shorter box
  double l00,l01,l10,l11;
  if(length(0,0) >= 0.0) l00 = length(0,0) + eps; else l00 = length(0,0) - eps;
  if(length(0,1) >= 0.0) l01 = length(0,1) + eps; else l01 = length(0,1) - eps;
  if(length(1,0) >= 0.0) l10 = length(1,0) + eps; else l10 = length(1,0) - eps;
  if(length(1,1) >= 0.0) l11 = length(1,1) + eps; else l11 = length(1,1) - eps;

  if (nlsd==2) 
    {
      // Rotate coordinates in tmp; Determine number of atoms
      for (int p=-2*ncells[1];p<2*ncells[1];p++) 
      	for (int q=-2*ncells[0];q<2*ncells[0];q++)    
	  {
	    dArrayT c(nlsd);
	    c[0] = (double)q; c[1] = (double)p;
	    for (int m=0;m<nuca;m++) 
	      {
		if ( natom >= temp_nat) {cout << "natoms wrong ";throw eSizeMismatch;}
		
		x = length(0,0);
		y = length(1,0);
		
		for (int k=0;k<nlsd;k++) 
		  {
		    x += (c[k] + vB(k,m))*vA(k,0);
		    y += (c[k] + vB(k,m))*vA(k,1);
		    type = vT[m];
		  }

		if(x >= l00 && x <= l01 && y >= l10 && y <= l11 )
		  {
		    (*temp_atom)(natom)[0] = x;
		    (*temp_atom)(natom)[1] = y;
		    (*temp_type)[natom] = type;
		    natom++;
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
		  if (natom > temp_nat) {cout << "natoms wrong ";throw eSizeMismatch;}

		  x = length(0,0);
		  y = length(1,0);
		  z = length(2,0);
		  
		  for (int k=0;k<nlsd;k++) 
		    {
		      x += (c[k] + vB(k,m))*vA(k,0);
		      y += (c[k] + vB(k,m))*vA(k,1);
		      z += (c[k] + vB(k,m))*vA(k,2);
		      type = vT[m];
		    }
		  
		  if(x >= l00 && x <= l01 && y >= l10 && y <= l11 &&
		     z >= l20 && z <= l21)
		    {
		      (*temp_atom)(natom)[0] = x;
		      (*temp_atom)(natom)[1] = y;
		      (*temp_atom)(natom)[2] = z;
		      (*temp_type)[natom] = type;
		      natom++;                     
		    }
		}
	    }
    }

  return natom;
}

int BoxT::RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat)
{
  int natom= 0;

  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();
  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vA = pcl->GetAxis();
  const dArray2DT& vB = pcl->GetBasis();
  const iArrayT& vT = pcl->GetType();

  if (nSD==2) 
    {
      for (int p=0;p<ncells[1];p++) 
	for (int q=0;q<ncells[0];q++) 
	  {
	    dArrayT c(nlsd);
	    c[0] = (double)q; c[1] = (double)p;

	    for (int m=0;m<nuca;m++) 
	      {
		if (natom > temp_nat) {cout << "natoms wrong";throw eSizeMismatch;}
		
		  (*temp_atom)(natom)[0] = length(0,0);
		  (*temp_atom)(natom)[1] = length(1,0);


		for (int k=0;k<nlsd;k++) 
		  {
		    (*temp_atom)(natom)[0] += (c[k] + vB(k,m))*vA(k,0);
		    (*temp_atom)(natom)[1] += (c[k] + vB(k,m))*vA(k,1);
		    (*temp_type)[natom] = vT[m];
		  }
		
		natom++;
	      }
	    
	  }
    }
  else if (nSD==3) 
    {
      for (int p=0;p<ncells[2];p++) 
	for (int q=0;q<ncells[1];q++) 
	  for (int r=0;r<ncells[0];r++) 
	    {
	      dArrayT c(nlsd);
	      c[0] = (double)r; c[1] = (double)q; c[2] = (double)p;
	      for (int m=0;m<nuca;m++) 
		{
		  
		  if (natom > temp_nat) {cout << "natoms wrong";throw eSizeMismatch;}
		  
		  (*temp_atom)(natom)[0] = length(0,0);
		  (*temp_atom)(natom)[1] = length(1,0);
		  (*temp_atom)(natom)[2] = length(2,0);
		  
		  for (int k=0;k<nlsd;k++) 
		    {
		      (*temp_atom)(natom)[0] += (c[k] + vB(k,m))*vA(k,0);
		      (*temp_atom)(natom)[1] += (c[k] + vB(k,m))*vA(k,1);
		      (*temp_atom)(natom)[2] += (c[k] + vB(k,m))*vA(k,2);
		      (*temp_type)[natom] = vT[m];
		    }
		  
		  natom++;        
		}
	      
	    }
    }

  return natom;
}


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

  for (int i=0; i < nSD; i++)
    {
      if (pbc[i]==1)
        {
          minmax(i,0) = length(i,0);
          minmax(i,1) = length(i,1);
        }
    }

  return minmax;  
}
