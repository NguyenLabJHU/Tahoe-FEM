/* $Id: BoxT.h,v 1.4 2002-08-02 02:07:49 saubry Exp $ */

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
  
 protected:
  
  iArrayT ncells;
  dArray2DT length; // lower and upper bounds
  
 public:
  
  //Constructor
  BoxT(int dim, dArray2DT len, dArrayT lattice_parameter);
  BoxT(int dim, iArrayT cel, dArrayT lattice_parameter);
  
  //Destructor
  ~BoxT(){};
  
  // Copy constructor
  BoxT(const BoxT& source);
  
  void CreateLattice(CrystalLatticeT* pcl); // return volume
  void CalculateBounds(iArrayT per,CrystalLatticeT* pcl);
  void CalculateType();

  iArrayT GetNCells();
  dArray2DT GetLength();

 private:

  dArray2DT ComputeMinMax();  

};

inline iArrayT BoxT::GetNCells(){return ncells;};
inline dArray2DT BoxT::GetLength(){return length;};

#endif
