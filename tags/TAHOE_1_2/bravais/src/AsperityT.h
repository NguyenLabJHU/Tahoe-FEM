// DEVELOPMENT
/* $Id: AsperityT.h,v 1.5 2003-06-13 17:27:28 saubry Exp $ */

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

 protected:

  double fRadius;
  dArrayT fCenterPlus;
  dArrayT fCenterMinus;

 public:
  
  //Constructor
  AsperityT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort);
  AsperityT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort);
  
  //Destructor
  ~AsperityT(){};
  
  // Copy constructor
  AsperityT(const AsperityT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); 
  void SortLattice(CrystalLatticeT* pcl);

  void CalculateBounds(iArrayT per,CrystalLatticeT* pcl);

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
