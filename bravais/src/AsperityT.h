// DEVELOPMENT
/* $Id: AsperityT.h,v 1.11 2003-09-08 20:14:51 jzimmer Exp $ */

#ifndef _ASPERITY_T_H_
#define _ASPERITY_T_H_

#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "VolumeT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class AsperityT : public VolumeT 
{
 protected:
 
  iArrayT ncells;
  dArray2DT length; // lower and upper bounds
  iArrayT WhichSort;
  iArrayT pbc;
  double fRadius;

 private:

  dArrayT fCenterPlus;
  dArrayT fCenterMinus;

  iArray2DT type1,type2;

 public:
  
  //Constructor
  AsperityT(int dim, dArray2DT len, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per);
  AsperityT(int dim, iArrayT cel, dArrayT lattice_parameter,
       iArrayT which_sort, StringT slt, iArrayT per);
  
  //Destructor
  ~AsperityT(){};
  
  // Copy constructor
  AsperityT(const AsperityT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); 
  void SortLattice(CrystalLatticeT* pcl);

  void CalculateBounds(CrystalLatticeT* pcl);

  iArrayT GetNCells();
  dArray2DT GetLength();
  double GetRadius();
 private: 

    dArray2DT ComputeMinMax(); 
    int RotateAtomInBox(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			iArrayT* temp_type, iArrayT* temp_parts,
			int temp_nat);
    int RotateBoxOfAtom(CrystalLatticeT* pcl,dArray2DT* temp_atom,
			iArrayT* temp_type, iArrayT* temp_parts, 
			int temp_nat);

    double ComputeCircleParameters();
};

inline double AsperityT::GetRadius() {return fRadius;};
inline iArrayT AsperityT::GetNCells(){return ncells;};
inline dArray2DT AsperityT::GetLength(){return length;};

#endif
