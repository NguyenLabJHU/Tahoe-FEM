// DEVELOPMENT
/* $Id: BoxT.h,v 1.21 2004-03-11 22:46:26 jzimmer Exp $ */

#ifndef _BOX_T_H_
#define _BOX_T_H_

#include <iostream.h>
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

 private:
  
  iArray2DT type1,type2;

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

    double Mod(double a,double p);
    int Maxx(int a,int b);
    double Maxx(double a,double b,double c);
    double Minn(double a,double b,double c);

    dArrayT CrossProduct(dArrayT,dArrayT);
    double DotProduct(dArrayT x,dArrayT y);

    double CalculatePeriodicLength(CrystalLatticeT* pcl,dArrayT Rot);

};

inline iArrayT BoxT::GetNCells(){return ncells;};
inline dArray2DT BoxT::GetLength(){return length;};

#endif
