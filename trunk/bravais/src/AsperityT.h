// DEVELOPMENT
/* $Id: AsperityT.h,v 1.6 2003-07-25 18:18:34 jzimmer Exp $ */

#ifndef _ASPERITY_T_H_
#define _ASPERITY_T_H_

#include <iostream>
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "VolumeT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class AsperityT : public VolumeT 
{
 
  iArrayT ncells;
  dArray2DT length; // lower and upper bounds
  iArrayT WhichSort;
  iArrayT pbc;

 protected:

  double fRadius;
  dArrayT fCenterPlus;
  dArrayT fCenterMinus;

 public:
  
  //Constructor
  AsperityT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort, iArrayT per);
  AsperityT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort, iArrayT per);
  
  //Destructor
  ~AsperityT(){};
  
  // Copy constructor
  AsperityT(const AsperityT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); 
  void SortLattice(CrystalLatticeT* pcl);

  void CalculateBounds(CrystalLatticeT* pcl);

  iArrayT GetNCells();
  dArray2DT GetLength();

 private: 

    dArray2DT ComputeMinMax(); 
    int RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			iArrayT* temp_parts, int temp_nat);
    int RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			iArrayT* temp_parts, int temp_nat);

    double ComputeCircleParameters();
};

inline iArrayT AsperityT::GetNCells(){return ncells;};
inline dArray2DT AsperityT::GetLength(){return length;};

#endif
