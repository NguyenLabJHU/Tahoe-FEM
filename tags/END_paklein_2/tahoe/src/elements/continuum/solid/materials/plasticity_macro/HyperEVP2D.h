/* $Id: HyperEVP2D.h,v 1.3.8.1 2002-10-28 06:49:24 paklein Exp $ */
#ifndef _HYPER_EVP_2D_H_
#define _HYPER_EVP_2D_H_

#include "HyperEVP3D.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class ElasticT;

class HyperEVP2D : public HyperEVP3D, public Material2DT
{
 public:
  // constructor
  HyperEVP2D(ifstreamT& in, const FDMatSupportT& support);

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
  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

} // namespace Tahoe 
#endif /* _HYPER_EVP_2D_ */
