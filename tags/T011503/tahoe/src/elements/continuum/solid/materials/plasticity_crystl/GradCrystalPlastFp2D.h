/* $Id: GradCrystalPlastFp2D.h,v 1.3 2002-11-14 17:06:32 paklein Exp $ */
#ifndef _GRAD_CRYSTAL_PLAST_FP_2D_H_
#define _GRAD_CRYSTAL_PLAST_FP_2D_H_

#include "GradCrystalPlastFp.h"
#include "Material2DT.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class GradCrystalPlastFp2D: public GradCrystalPlastFp, public Material2DT
{
 public:
  // constructor
  GradCrystalPlastFp2D(ifstreamT& in, const FDMatSupportT& support);

  // destructor
  ~GradCrystalPlastFp2D();

  // crystal Cauchy stress
  virtual const dSymMatrixT& s_ij();

  // crystal modulus 
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected: 

  // crystal Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dc_ijkl;
};

} // namespace Tahoe 
#endif /* _GRAD_CRYSTAL_PLAST_FP_2D_H_ */
