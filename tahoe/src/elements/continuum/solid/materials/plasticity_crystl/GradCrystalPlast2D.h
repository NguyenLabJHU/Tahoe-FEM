/*
  File: GradCrystalPlast2D.h
*/

#ifndef _GRAD_CRYSTAL_PLAST_2D_H_
#define _GRAD_CRYSTAL_PLAST_2D_H_

#include "GradCrystalPlast.h"
#include "Material2DT.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

class GradCrystalPlast2D: public GradCrystalPlast, public Material2DT
{
 public:
  // constructor
  GradCrystalPlast2D(ifstreamT& in, const ElasticT& element);

  // destructor
  ~GradCrystalPlast2D();

  // crystal Cauchy stress
  virtual const dSymMatrixT& s_ij();

  // crystal modulus 
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected: 

//DEV - deprecated
  // deformation gradient 
  //virtual const dMatrixT& DeformationGradient(const LocalArrayT& disp);

 protected: 
  // crystal Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dc_ijkl;

};

#endif /* _GRAD_CRYSTAL_PLAST_2D_H_ */
