/* $Id: GradCrystalPlastFp2D.h,v 1.6 2004-09-10 22:39:43 paklein Exp $ */
#ifndef _GRAD_CRYSTAL_PLAST_FP_2D_H_
#define _GRAD_CRYSTAL_PLAST_FP_2D_H_

#include "GradCrystalPlastFp.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class GradCrystalPlastFp2D: public GradCrystalPlastFp
{
 public:
  // constructor
  GradCrystalPlastFp2D(ifstreamT& in, const FSMatSupportT& support);

  // crystal Cauchy stress
  virtual const dSymMatrixT& s_ij();

  // crystal modulus 
  virtual const dMatrixT& c_ijkl();

 protected: 

  // crystal Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dc_ijkl;
};

} // namespace Tahoe 
#endif /* _GRAD_CRYSTAL_PLAST_FP_2D_H_ */
