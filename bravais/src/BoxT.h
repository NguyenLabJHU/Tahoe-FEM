// DEVELOPMENT
/* $Id: BoxT.h,v 1.13 2003-07-21 15:27:58 fwdelri Exp $ */

#ifndef _BOX_T_H_
#define _BOX_T_H_

#include <iostream>
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "VolumeT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class BoxT : public VolumeT 
{
  
  iArrayT ncells;
  dArray2DT length; // lower and upper bounds
  iArrayT WhichSort;

 public:
  
  //Constructor
  BoxT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt);
  BoxT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt);
  
  //Destructor
  ~BoxT(){};
  
  // Copy constructor
  BoxT(const BoxT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); 
  void SortLattice(CrystalLatticeT* pcl);

  void CalculateBounds(iArrayT per,CrystalLatticeT* pcl);

  iArrayT GetNCells();
  dArray2DT GetLength();

 private:

    dArray2DT ComputeMinMax(); 
    int RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat);
    int RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,iArrayT* temp_type,int temp_nat);

};

inline iArrayT BoxT::GetNCells(){return ncells;};
inline dArray2DT BoxT::GetLength(){return length;};

#endif
