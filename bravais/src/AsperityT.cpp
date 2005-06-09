// DEVELOPMENT
/* $Id: AsperityT.cpp,v 1.19 2005-06-09 21:48:08 saubry Exp $ */
#include "AsperityT.h"
#include "VolumeT.h"

#include "ExceptionCodes.h"
#include "ifstreamT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

AsperityT::AsperityT(int dim, dArray2DT len,
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
  VolType = "ASPERITY";

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
	  /*
          int ncl = static_cast<int>(len(i,0)/lattice_parameter[i]);
	  int ncu = static_cast<int>(len(i,1)/lattice_parameter[i]);
          length(i,0) = ncl * lattice_parameter[i];
          length(i,1) = ncu * lattice_parameter[i];

	  ncells[i] = ncu - 1 - ncl;
	  */

	  double dist = ncells[i]*lattice_parameter[i]*0.5;       
	  length(i,0) = -dist;       
	  length(i,1) = length(i,0) + (dist - length(i,0));
        }
    }
}

AsperityT::AsperityT(int dim, iArrayT cel,
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
  VolType = "ASPERITY";

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

AsperityT::AsperityT(const AsperityT& source) : VolumeT(source.nSD) 
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

  atom_coord.Dimension(source.nATOMS,source.nSD);
  atom_coord = source.atom_coord;

  atom_number.Dimension(source.nATOMS);
  atom_number = source.atom_number;

  atom_types.Dimension(source.nATOMS);
  atom_types = source.atom_types;

  atom_ID.Dimension(source.atom_ID.Length());
  atom_ID = source.atom_ID;

  atom_connect.Dimension(source.atom_connect.Length());
  atom_connect = source.atom_connect;

  VolType = "ASPERITY";
}

void AsperityT::CreateLattice(CrystalLatticeT* pcl) 
{
  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();
  int ntype = pcl->GetNTYPE();

  int temp_nat=0;
  dArray2DT temp_atom;
  iArrayT temp_type;
  iArrayT temp_parts;

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
  temp_parts.Dimension(temp_nat);
  temp_type.Dimension(temp_nat);

  if(pcl->GetRotMeth() == 0)
    nATOMS = RotateAtomInBox(pcl,&temp_atom,&temp_type,&temp_parts,temp_nat);
  else
    nATOMS = RotateBoxOfAtom(pcl,&temp_atom,&temp_type,&temp_parts,temp_nat);

  // Get atoms coordinates
  atom_coord.Dimension(nATOMS,nlsd);
  atom_number.Dimension(nATOMS);
  atom_types.Dimension(nATOMS);

  for(int m=0; m < nATOMS ; m++)
    {
      atom_number[m] = m;
      atom_types[m] = temp_type[m];
      for (int k=0;k< nlsd;k++)
	atom_coord(m)[k] = temp_atom(m)[k];
    }

  atom_names = "Asperity";

  // Create parts
  atom_parts.Dimension(nATOMS);
  for(int m=0; m < nATOMS ; m++) 
    atom_parts[m] = temp_parts[m];

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

  // Create types and connectivities
  if(ntype > 2) 
    cout << "WARNING: ** nTypes == 2  maximum for ensight output ** \n";

  atom_ID.Dimension(ntype);
  atom_connect.Dimension(ntype);

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

  atom_ID[0] = "type1";
  if (ntype > 1) atom_ID[1] = "type2";

  if (n < nATOMS && n > 0) type1.Resize(n);
  if (m < nATOMS && m > 0) type2.Resize(m);
  atom_connect[0] = &type1;
  if (ntype > 1) atom_connect[1] = &type2;

}

void AsperityT::SortLattice(CrystalLatticeT* pcl) 
{
  int nlsd = pcl->GetNLSD();
  
  dArray2DT new_coord;
  new_coord.Dimension(atom_coord.MajorDim(),nlsd);   
  new_coord = 0.0;

  iArrayT new_type;
  new_type.Dimension(nATOMS);  
  new_type = 0;

  iArrayT typ(atom_coord.MajorDim());typ = 0;

  dArrayT x(atom_coord.MajorDim()); x = 0.0;
  dArrayT y(atom_coord.MajorDim()); y = 0.0;
  dArrayT z(atom_coord.MajorDim()); z = 0.0;

  iArrayT part(atom_parts.Length()); part = 0;

  iArrayT Map(atom_coord.MajorDim()); Map = 0;

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
      part[m] = atom_parts[Map[m]];
    } 

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;
  atom_parts = part;

  // Sort 2nd criterium
  new_coord = 0.0;
  new_type = 0;
  part = 0;

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

  x = 0.0;
  y = 0.0;
  z = 0.0;
  typ = 0;
  int p = 0;
  for(int n = 0; n < is-1; n++)
    {
      dArrayT aux(Ind[n+1]-Ind[n]);
      dArrayT aux2(Ind[n+1]-Ind[n]);
      iArrayT aux3(Ind[n+1]-Ind[n]);
      aux = 0.0;
      aux2= 0.0;
      iArrayT aux_par(Ind[n+1]-Ind[n]); aux_par = 0;

      int isa = 0;
      for(int m = Ind[n]; m < Ind[n+1]; m++)
	{
	  aux[isa] = atom_coord(m)[WhichSort[1]];
	  if(nlsd == 3) aux2[isa]= atom_coord(m)[WhichSort[2]];
	  aux3[isa] = atom_types[m];
	  aux_par[isa] = atom_parts[m];
	  isa++;
	}

      iArrayT Map2(isa);
      Map2.SetValueToPosition();
      Map2.SortAscending(aux);

      for(int m = 0; m < isa; m++)
	{
	  y[p] = aux[m];
	  if(nlsd == 3) z[p] = aux2[Map2[m]];
	  typ[p] = aux3[Map2[m]];
	  part[p] = aux_par[Map2[m]];
	  p++;
	}
    }
  
  for(int m=0; m < nATOMS ; m++) 
    {
      new_coord(m)[WhichSort[0]] = atom_coord(m)[WhichSort[0]];
      new_coord(m)[WhichSort[1]] = y[m];
      if (nlsd == 3) new_coord(m)[WhichSort[2]] = z[m];
      new_type[m] = typ[m];
    } 

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;
  atom_parts = part;

  // Sort 3nd criterium
  if(nlsd == 3)
  {
    new_coord = 0.0;
    new_type = 0;
    part = 0;

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
    
    x = 0.0;
    y = 0.0;
    z = 0.0;
    typ = 0;
    int p = 0;
    for(int n = 0; n < is-1; n++)
      {
	dArrayT aux(Ind[n+1]-Ind[n]);
	aux = 0.0;
	iArrayT aux2(Ind[n+1]-Ind[n]);
	aux2 = 0;
	iArrayT aux_par(Ind[n+1]-Ind[n]); aux_par = 0;

	int isa = 0;
	for(int m = Ind[n]; m < Ind[n+1]; m++)
	  {
	    aux[isa] = atom_coord(m)[WhichSort[2]];
	    aux2[isa] = atom_types[m];
	    aux_par[isa] = atom_parts[m];
	    isa++;
	  }
	
	iArrayT Map2(isa);
	Map2.SetValueToPosition();
	Map2.SortAscending(aux);
	
	for(int m = 0; m < isa; m++)
	  {
	    z[p] = aux[m];
	    typ[p] = aux2[m];
	    part[p] = aux_par[Map2[m]];
	    p++;
	  }
      }
    
    for(int m=0; m < nATOMS ; m++) 
      {
	new_coord(m)[WhichSort[0]] = atom_coord(m)[WhichSort[0]];
	new_coord(m)[WhichSort[1]] = atom_coord(m)[WhichSort[1]];
	new_coord(m)[WhichSort[2]] = z[m];
	new_type[m] = typ[m];
      } 
  }

  // Update sorted atoms
  atom_coord = new_coord;
  atom_types = new_type;
  atom_parts = part;
}

void AsperityT::CalculateBounds(CrystalLatticeT* pcl)
{
  const dArrayT& vLP = pcl->GetLatticeParameters();
  atom_bounds.Dimension(nSD,2);

  for (int i=0; i < nSD; i++)
    {
      if (pbc[i]==0) 
	{
	  // non-periodic conditions
          atom_bounds(i,0) = length(i,0);
          atom_bounds(i,1) = length(i,1);
	  //atom_bounds(i,0) = -10000.; -> different for Alex's code.
	  //atom_bounds(i,1) =  10000.;
	}
      else if (pbc[i]==1)
	{
	  // periodic conditions
          atom_bounds(i,0) = length(i,0) + 0.75 * vLP[i]; 
	  // takes the shift into account.
          atom_bounds(i,1) = length(i,1) + 0.75 * vLP[i];
	}
      else
	throw eBadInputValue;
    }  
}


//////////////////// PRIVATE //////////////////////////////////

double AsperityT::ComputeCircleParameters()
{
  double rx,rz,h0;

  fCenterPlus.Dimension(nSD);  
  fCenterMinus.Dimension(nSD);  

  h0 = 0.0;  //length(nSD-1,0) + (length(nSD-1,1) - length(nSD-1,0))*0.5;
  rx = (length(0,1)-length(0,0))*0.5;
  rz = length(nSD-1,1)-h0;
  fRadius = (rx*rx+rz*rz)*0.5/rz;

  fCenterPlus[0] = length(0,0) + (length(0,1) - length(0,0))*0.5;
  if(nSD == 3)
    fCenterPlus[1] = length(1,0) + (length(1,1) - length(1,0))*0.5;
  fCenterPlus[nSD-1] = h0 + fRadius;

  fCenterMinus[0] = length(0,0) + (length(0,1) - length(0,0))*0.5;
  if(nSD == 3)
    fCenterMinus[1] = length(1,0) + (length(1,1) - length(1,0))*0.5;
  fCenterMinus[nSD-1] = h0 - fRadius;

  return h0;
}

int AsperityT::RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			       iArrayT* temp_type,iArrayT* temp_parts,int temp_nat)
{
  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();

  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vB = pcl->GetBasis();
  const dArray2DT& vA = pcl->GetAxis();
  const iArrayT& vT = pcl->GetType();

  double x,y,z=0;
  double eps = 1.e-6;

  // Call circle parameters
  double h0 = ComputeCircleParameters();
  cout << "(RAB) Height of the asperity is " << h0 << "\n";
  cout << "(RAB) Radius of the asperity is " << fRadius << "\n";

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
		    
		double r0 = x - fCenterPlus[0]; 
		double r1 = y - fCenterPlus[1]; 
		double R = r0*r0 + r1*r1;

		if(x >= l00 && x <= l01 && y >= l10 && y <= l11)
		   if(R <= fRadius*fRadius || z <= h0 )
		  {
		    (*temp_atom)(natom)[0] = x;
		    (*temp_atom)(natom)[1] = y;
		    (*temp_type)[natom] = type;
		    (*temp_parts)[natom] = 1; 
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

		  // Asperity against a block
		  double r0 = x - fCenterPlus[0]; 
		  double r1 = y - fCenterPlus[1]; 
		  double r2 = z - fCenterPlus[2]; 
		  double R = r0*r0 + r1*r1 + r2*r2;

		  if(x >= l00 && x <= l01 && y >= l10 && y <= l11 && z >= l20 && z <= l21)
		    if(R <= fRadius*fRadius || z < h0)
		      {
			if( fabs(x-length(0,0)) >= 1.e-5 && fabs(y-length(1,0)) >= 1.e-5 ) 
			  {
			    (*temp_atom)(natom)[0] = x;
			    (*temp_atom)(natom)[1] = y;
			    (*temp_atom)(natom)[2] = z;
			    (*temp_type)[natom] = type;
			    (*temp_parts)[natom] = 1; 
			    if (z< h0) (*temp_parts)[natom] = -1; 
			    if ( fabs(z-length(2,0)) <= 1.e-5 ) 
			      (*temp_parts)[natom]= -2;
			    
			    natom++;                     
			  }
		      }
		 

		  // Two Asperities
		  /*
		  double rp0 = x - fCenterPlus[0]; 
		  double rp1 = y - fCenterPlus[1]; 
		  double rp2 = z - fCenterPlus[2]; 
		  double Rplus = rp0*rp0 + rp1*rp1 + rp2*rp2;

		  double rm0 = x - fCenterMinus[0]; 
		  double rm1 = y - fCenterMinus[1]; 
		  double rm2 = z - fCenterMinus[2]; 
		  double Rmin = rm0*rm0 + rm1*rm1 + rm2*rm2;

		  if(x >= l00 && x <= l01 && y >= l10 && y <= l11 && z >= l20 && z <= l21)
		    if( (Rplus <= fRadius*fRadius && z >= h0) || 
		        (Rmin  <= fRadius*fRadius && z <= h0) )
		      {
			(*temp_atom)(natom)[0] = x;
			(*temp_atom)(natom)[1] = y;
			(*temp_atom)(natom)[2] = z;
			(*temp_parts)[natom] = 10;
			if (Rmin <= fRadius*fRadius && z <= h0 + eps) 
			     (*temp_parts)[natom] = 20;
			if ( fabs(z-length(2,0)) <= 1.e-5 ) 
			     (*temp_parts)[natom]= 30;
			natom++;                     
		      }
		  */
		}
		  
	    }


      // Shift asperity
      for (int k=0;k<natom;k++) 
	{
	  if((*temp_parts)[k] > 0) 
	    {
	      (*temp_atom)(k)[0] += vLP[0]*0.25;
	      (*temp_atom)(k)[1] += vLP[1]*0.25;
	    }
	}


    }
  return natom;
}

int AsperityT::RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			       iArrayT* temp_type,iArrayT* temp_parts,int temp_nat)
{
  int natom= 0;

  int nlsd = pcl->GetNLSD();
  int nuca = pcl->GetNUCA();
  const dArrayT& vLP = pcl->GetLatticeParameters();
  const dArray2DT& vA = pcl->GetAxis();
  const dArray2DT& vB = pcl->GetBasis();
  const iArrayT& vT = pcl->GetType();

  double eps = 1.e-6;
  double x,y,z=0;

  // Call circle parameters
  double h0 = ComputeCircleParameters();
  cout << "(RBA) Height of the asperity is " << h0 << "\n";
  cout << "(RBA) Radius of the asperity is " << fRadius << "\n";

  dArrayT rotated_fCenterPlus(nlsd);
  rotated_fCenterPlus = pcl->VectorRotation(fCenterPlus);

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
		
		x = length(0,0);
		y = length(1,0);

		for (int k=0;k<nlsd;k++) 
		  {
		    x += (c[k] + vB(k,m))*vA(k,0);
		    y += (c[k] + vB(k,m))*vA(k,1);
		  }

		double r0 = x - rotated_fCenterPlus[0]; 
		double r1 = y - rotated_fCenterPlus[1]; 
		double R = r0*r0 + r1*r1;

		if(R <= fRadius*fRadius || z <= h0 )
		  {
		    (*temp_atom)(natom)[0] = x;
		    (*temp_atom)(natom)[1] = y;
		    (*temp_type)[natom] = vT[m];
		    (*temp_parts)[natom] = 1; 
		    natom++;
		  }
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
		  
		  x = length(0,0);
		  y = length(1,0);
		  z = length(2,0);
		  
		  for (int k=0;k<nlsd;k++) 
		    {
		      x += (c[k] + vB(k,m))*vA(k,0);
		      y += (c[k] + vB(k,m))*vA(k,1);
		      z += (c[k] + vB(k,m))*vA(k,2);
		    }
		  
		  // Asperity against a block
		  double r0 = x - rotated_fCenterPlus[0]; 
		  double r1 = y - rotated_fCenterPlus[1]; 
		  double r2 = z - rotated_fCenterPlus[2]; 
		  double R = r0*r0 + r1*r1 + r2*r2;
		  
		  if( (R + eps < fRadius*fRadius ) || (z < h0) ) 
		    {
		      (*temp_atom)(natom)[0] = x;
		      (*temp_atom)(natom)[1] = y;
		      (*temp_atom)(natom)[2] = z;
		      (*temp_type)[natom] = vT[m];
		      (*temp_parts)[natom] = 1;
		      if (z< h0) (*temp_parts)[natom] = -1; 
		      if ( fabs(z-length(2,0)) <= 1.e-5 ) 
			(*temp_parts)[natom]= -2;
		      natom++;        
		    }

		  // Two Asperities
		  /*
		  double rp0 = x - fCenterPlus[0]; 
		  double rp1 = y - fCenterPlus[1]; 
		  double rp2 = z - fCenterPlus[2]; 
		  double Rplus = rp0*rp0 + rp1*rp1 + rp2*rp2;

		  double rm0 = x - fCenterMinus[0]; 
		  double rm1 = y - fCenterMinus[1]; 
		  double rm2 = z - fCenterMinus[2]; 
		  double Rmin = rm0*rm0 + rm1*rm1 + rm2*rm2;

		  if(x >= l00 && x <= l01 && y >= l10 && y <= l11 && z >= l20 && z <= l21)
		    if( (Rplus <= fRadius*fRadius && z >= h0) || 
		        (Rmin  <= fRadius*fRadius && z <= h0) )
		      {
			(*temp_atom)(natom)[0] = x;
			(*temp_atom)(natom)[1] = y;
			(*temp_atom)(natom)[2] = z;
			(*temp_type)[natom] = vT[m];
			(*temp_parts)[natom] = 1;
			if (Rmin <= fRadius*fRadius && z <= h0 + eps) 
			     (*temp_parts)[natom] = -1;
			if ( fabs(z-length(2,0)) <= 1.e-5 ) 
			     (*temp_parts)[natom]= 2;
			natom++;                     
		      }
		  */

		}
	      
	    }

      // Shift asperity
      for (int k=0;k<natom;k++) 
	{
	  if((*temp_parts)[k] > 0) 
	    {
	      (*temp_atom)(k)[0] += vLP[0]*0.25;
	      (*temp_atom)(k)[1] += vLP[1]*0.25;
	    }
	}

    }
  return natom;
}


dArray2DT AsperityT::ComputeMinMax()
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
      if (pbc[i]==1)
        {
          minmax(i,0) = length(i,0);
          minmax(i,1) = length(i,1);
        }
    }


  /*  for (int i=0; i < nSD; i++)
    {
      if(minmax(i,0) > minmax(i,1)) 
	{
	  double temp = minmax(i,0);
	  minmax(i,0) = minmax(i,1);
	  minmax(i,1) = temp;
	}
    }
  */

  return minmax;  
}
