/*
  File: BCJHypoIsoDamageYC2D.h
*/

#ifndef _BCJ_HYPO_ISO_DAMAGE_YC_2D_H_
#define _BCJ_HYPO_ISO_DAMAGE_YC_2D_H_

#include "BCJHypoIsoDamageYC3D.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

class ifstreamT;
class ElasticT;

class BCJHypoIsoDamageYC2D : public BCJHypoIsoDamageYC3D, public Material2DT
{
 public:
  // constructor
  BCJHypoIsoDamageYC2D(ifstreamT& in, const FiniteStrainT& element);

  // destructor
  ~BCJHypoIsoDamageYC2D();

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:

  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

#endif /* _BCJ_HYPO_ISO_DAMAGE_YC_2D_ */
