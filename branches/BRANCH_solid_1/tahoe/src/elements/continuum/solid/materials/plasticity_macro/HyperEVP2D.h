/*
  File: HyperEVP2D.h
*/

#ifndef _HYPER_EVP_2D_H_
#define _HYPER_EVP_2D_H_

#include "HyperEVP3D.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

class ifstreamT;
class ElasticT;

class HyperEVP2D : public HyperEVP3D, public Material2DT
{
 public:
  // constructor
  HyperEVP2D(ifstreamT& in, const ElasticT& element);

  // destructor
  ~HyperEVP2D();

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:
  // deformation gradient
  //virtual const dMatrixT& DeformationGradient(const LocalArrayT& disp);
//DEV - see who needs this

 protected:
  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

#endif /* _HYPER_EVP_2D_ */
