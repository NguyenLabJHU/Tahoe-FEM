// DEVELOPMENT
/* $Id: BoxT.h,v 1.14 2003-07-25 18:18:34 jzimmer Exp $ */

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
  iArrayT pbc;

 public:
  
  //Constructor
  BoxT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per);
  BoxT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per);
  
  //Destructor
  ~BoxT(){};
  
  // Copy constructor
  BoxT(const BoxT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); 
  void SortLattice(CrystalLatticeT* pcl);

  void CalculateBounds(CrystalLatticeT* pcl);

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
