// DEVELOPMENT
/* $Id: BoxT.cpp,v 1.47 2005-02-11 18:33:11 saubry Exp $ */
#include "BoxT.h"
#include "VolumeT.h"

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
      VolType = "BOX";
   
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
	   //double dist = ncells[i]*lattice_parameter[i]*0.5;
            //length(i,0) = 0 ; // -dist;
	    // length(i,1) = dist; // length(i,0) + (dist - length(i,0));

	    length(i,0) = len(i,0);
	    length(i,1) = len(i,1);

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
      VolType = "BOX";
   
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
   
      atom_coord.Dimension(source.nATOMS,source.nSD);
      atom_coord = source.atom_coord;
   
      atom_types.Dimension(source.nATOMS);
      atom_types = source.atom_types;
   
      atom_ID.Dimension(source.atom_ID.Length());
      atom_ID = source.atom_ID;
   
      atom_connect.Dimension(source.atom_connect.Length());
      atom_connect = source.atom_connect;
   
      VolType = "BOX";
   }


    void BoxT::CreateLattice(CrystalLatticeT* pcl) 
   {
      int nlsd = pcl->GetNLSD();
      int nuca = pcl->GetNUCA();
      int ntype = pcl->GetNTYPE();
   
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
   
      atom_names = "Box";
   
   // Update lengths
      length = ComputeMinMax();
   
   //Calculate volume here
      switch(nlsd) 
	{
	case 2:
	  volume = (length(0,1)-length(0,0))*(length(1,1)-length(1,0));
	  break;
	case 3:
	  volume = (length(0,1)-length(0,0))*         
	    (length(1,1)-length(1,0))*
               (length(2,1)-length(2,0));
	  break;
	}

   // Output spacings along ox, oy and oz.

       iArrayT TmpSort;
       TmpSort.Dimension(nSD);
       TmpSort = WhichSort;

       WhichSort[0] = 1;
       WhichSort[1] = 2;
       WhichSort[2] = 0;
       SortLattice(pcl);

       WhichSort[0] = 2;
       WhichSort[1] = 0;
       WhichSort[2] = 1;
       SortLattice(pcl);

       WhichSort[0] = 0;
       WhichSort[1] = 1;
       WhichSort[2] = 2;
       SortLattice(pcl);
       cout << "\n";

       WhichSort = TmpSort;
       
   
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
   
   
   // Create parts
      atom_parts.Dimension(nATOMS);
      atom_parts = 1;
   
   }

    void BoxT::SortLattice(CrystalLatticeT* pcl)
   {
      int nlsd = pcl->GetNLSD();
   
      dArray2DT new_coord;
      new_coord.Dimension(atom_coord.MajorDim(),nlsd);   
      new_coord = 0.0;
   
      dArrayT x(atom_coord.MajorDim()); x = 0.0;
      dArrayT y(atom_coord.MajorDim()); y = 0.0;
      dArrayT z(atom_coord.MajorDim()); z = 0.0;
   
      iArrayT new_type;
      new_type.Dimension(atom_coord.MajorDim());  
      new_type = 0;
   
      iArrayT typ(atom_coord.MajorDim());typ = 0;
   
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
      
         int isa = 0;
         for(int m = Ind[n]; m < Ind[n+1]; m++)
         {
            aux[isa] = atom_coord(m)[WhichSort[1]];
            if(nlsd == 3) aux2[isa]= atom_coord(m)[WhichSort[2]];
            aux3[isa] = atom_types[m];
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
      
         x = 0.0;
         y = 0.0;
         z = 0.0;
         typ = 0;
         int p = 0;
         for(int n = 0; n < is-1; n++)
         {
            dArrayT aux(Ind[n+1]-Ind[n]);
            iArrayT aux2(Ind[n+1]-Ind[n]);
            aux = 0.0;
            aux2 = 0;
         
            int isa = 0;
            for(int m = Ind[n]; m < Ind[n+1]; m++)
            {
               aux[isa] = atom_coord(m)[WhichSort[2]];
               aux2[isa] = atom_types[m];
               isa++;
            }
         
            iArrayT Map2(isa);
            Map2.SetValueToPosition();
            Map2.SortAscending(aux);
         
            for(int m = 0; m < isa; m++)
            {
               z[p] = aux[m];
               typ[p] = aux2[m];
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

      
      cout << "Periodic distance in " << WhichSort[2] << ": " << z[1] - z[0] << "\n";

   
   // Update sorted atoms
      atom_coord = new_coord;
      atom_types = new_type;   
   }

    void BoxT::CalculateBounds(CrystalLatticeT* pcl)
   {
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
            atom_bounds(i,0) = length(i,0);
            atom_bounds(i,1) = length(i,1);

	    cout << "Bounds:";
	    cout << i << " " << atom_bounds(i,0) << " " <<  atom_bounds(i,1) << "\n";

         }
         else
            throw eBadInputValue;
      } 

 
   }


//////////////////// PRIVATE //////////////////////////////////

    int BoxT::RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat)
   {
      const dArrayT& vLP = pcl->GetLatticeParameters();

      int nlsd = pcl->GetNLSD();
      const dArrayT& vec_a = pcl->GetVector_a();
      const dArrayT& vec_b = pcl->GetVector_b();
      const dArrayT& vec_c = pcl->GetVector_c();
   
      double x,y,z;
         
      int natom= 0;
      int type= 0;
   
   // Define new number of cells after rotation
      double xlen = Maxx(fabs(vec_a[0]),fabs(vec_a[1]),fabs(vec_a[2]))*vLP[0]/2.;
      double ylen = Maxx(fabs(vec_b[0]),fabs(vec_b[1]),fabs(vec_b[2]))*vLP[1]/2.;
      double zlen = Maxx(fabs(vec_c[0]),fabs(vec_c[1]),fabs(vec_c[2]))*vLP[2]/2.;
   
      iArrayT new_cel;
      new_cel.Dimension(nlsd);
      new_cel[0] = (int)((length(0,1)-length(0,0))/xlen) + 1;
      new_cel[1] = (int)((length(1,1)-length(1,0))/ylen) + 1;
      new_cel[0] = Maxx(new_cel[0],new_cel[1]);
      if(nlsd == 3) 
      {
         new_cel[2] = (int)((length(2,1)-length(2,0))/zlen) + 1;
         new_cel[0] = Maxx(new_cel[0],new_cel[2]);
      }
      new_cel[0] = 2*new_cel[0];
      int ncentx = new_cel[0]/2;
   
      double epsx = 0.01*(length(0,1)-length(0,0));
      double epsy = 0.01*(length(1,1)-length(1,0));
      double l00 = length(0,0);//+ epsx;
      double l01 = length(0,1)- epsx;
      
      double l10 = length(1,0);//+ epsy;
      double l11 = length(1,1)- epsy;

      if (nlsd==2) 
      {
	cout << epsx << " " << epsy << "\n";
	// Rotate coordinates
	for (int r=0;r<new_cel[0];r++) 
	  for (int q=0;q<new_cel[0];q++)    	
	    {
	      if (natom > temp_nat) {cout << "natoms wrong ";throw eSizeMismatch;}
	      
	      int rr = r - ncentx;
	      int qq = q - ncentx;
               
	      x = (rr*vec_a[0] + qq*vec_b[0])*vLP[0]/2.;
	      y = (rr*vec_a[1] + qq*vec_b[1])*vLP[1]/2.;
	      
	      type = 1;
	      
	      if(x > l00 && x < l01 && y > l10 && y < l11)
		{
		  (*temp_atom)(natom)[0] = x;
		  (*temp_atom)(natom)[1] = y;
		  (*temp_type)[natom] = type;
		  natom++;                     
		}
	    }
	
      }
      else if (nlsd==3) 
      {
         double epsz = 0.01*(length(2,1)-length(2,0));      
         double l20 = length(2,0);//+ epsz;
         double l21 = length(2,1)- epsz;

	 cout << epsx << " " << epsy << " " << epsz << "\n";

      
      // Rotate coordinates
         for (int r=0;r<new_cel[0];r++) 
            for (int q=0;q<new_cel[0];q++)    
               for (int p=0;p<new_cel[0];p++) 	
               {
                  if (natom > temp_nat) {cout << "natoms wrong ";throw eSizeMismatch;}
          
                  int rr = r - ncentx;
                  int qq = q - ncentx;
                  int pp = p - ncentx;
               
                  x = (rr*vec_a[0] + qq*vec_b[0] + pp*vec_c[0])*vLP[0]/2.;
                  y = (rr*vec_a[1] + qq*vec_b[1] + pp*vec_c[1])*vLP[1]/2.;
                  z = (rr*vec_a[2] + qq*vec_b[2] + pp*vec_c[2])*vLP[2]/2.;
               
                  type = 1;
               
                  if(x > l00 && x < l01 && y > l10 && y < l11 &&
                  z > l20 && z < l21)
                  {
                     (*temp_atom)(natom)[0] = x;
                     (*temp_atom)(natom)[1] = y;
                     (*temp_atom)(natom)[2] = z;
                     (*temp_type)[natom] = type;
                     natom++;                     
                  }
               }
      }
      cout << "natom=" << natom <<"\n";
      return natom;
   }

    int BoxT::RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat)
   {
      int natom= 0;
   
      int nlsd = pcl->GetNLSD();
      int nuca = pcl->GetNUCA();
      const dArray2DT& vA = pcl->GetAxis();
      const dArray2DT& vB = pcl->GetBasis();
      const iArrayT& vT = pcl->GetType();
   
   
   //cout << "Rotating Box of Atoms\n";
   
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

    double BoxT::CalculatePeriodicLength(CrystalLatticeT* pcl,dArrayT Rot)
   {
      double PerLen = 0.0;
      double tol = 1.e-5;
      int nlsd = pcl->GetNLSD();
   
      const dArrayT& vec_a = pcl->GetVector_a();
      const dArrayT& vec_b = pcl->GetVector_b();
      const dArrayT& vec_c = pcl->GetVector_c();
      const dArrayT& vLP = pcl->GetLatticeParameters();
   
      dArrayT p(nlsd);
      double rmax = 100.0;
   
      dArrayT capa(nlsd);
      capa = CrossProduct(vec_b,vec_c);
      dArrayT capb(nlsd);
      capb = CrossProduct(vec_c,vec_a);
      dArrayT capc(nlsd);
      capc = CrossProduct(vec_a,vec_b);
   
      double rl = DotProduct(capa,Rot);
      double rm = DotProduct(capb,Rot);
      double rn = DotProduct(capc,Rot);
   
      double ql = fabs(rl);
      double qm = fabs(rm);
      double qn = fabs(rn);
   
   // Find smallest non-zero
      if(ql<tol) ql = 1.e6;
      if(qm<tol) qm = 1.e6;
      if(qn<tol) qn = 1.e6;
      double small = Minn(ql,qm,qn);
   
   // Multiply by 1/small
      rl /= small;
      rm /= small;
      rn /= small;
   
      ql=rl;
      qm=rm;
      qn=rn;
   
      double rmult = 1.0;
   // Check to see if integer
      while (rmult < rmax)
      {
         int Integer = 1;
         double tql=fabs(ql);
         double tqm=fabs(qm);
         double tqn=fabs(qn);
      
         if(fabs(Mod(tql+0.5,1.0)-0.5) > tol) Integer = 0;
         if(fabs(Mod(tqm+0.5,1.0)-0.5) > tol) Integer = 0;
         if(fabs(Mod(tqn+0.5,1.0)-0.5) > tol) Integer = 0;
      
         if(Integer == 1) 
         {
            for (int j=0; j < nlsd; j++)
               p[j] = ql*vec_a[j]*vLP[0]/2. + qm*vec_b[j]*vLP[1]/2. + qn*vec_c[j]*vLP[2]/2.;
            PerLen = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
            return PerLen;
         }
         else
         {
            ql=rmult*rl;
            qm=rmult*rm;
            qn=rmult*rn;
            rmult = rmult + 1.0;
         }
      }
   
      for (int j=0; j < nlsd; j++)
         p[j] = rl*vec_a[j] + rm*vec_b[j] + rn*vec_c[j];
      PerLen = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
      return PerLen;
   }

    double BoxT::Mod(double a,double p)
   {
      double b;
      b = a - (int)(a/p)*p;
      return b;
   }


    dArrayT BoxT::CrossProduct(dArrayT x,dArrayT y)
   {
      dArrayT z(x.Length());
   
      z[0] = x[1]*y[2] - x[2]*y[1];
      z[1] = x[2]*y[0] - x[0]*y[2];
      z[2] = x[0]*y[1] - x[1]*y[0];
      return z;
   }

    double BoxT::DotProduct(dArrayT x,dArrayT y)
   {
      return  x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
   }

    int BoxT::Maxx(int a,int b)
   {
      return (a>b)?a:b;
   }

    double BoxT::Maxx(double a,double b,double c)
   {
      double maxi;
      maxi = (a>b)?a:b;
      if(maxi<c) maxi = c;
      return maxi;
   }

    double BoxT::Minn(double a,double b,double c)
   {
      double mini;
      mini = (a<b)?a:b;
      if(mini>c) mini = c;
      return mini;
   }
