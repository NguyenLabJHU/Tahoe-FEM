/* $Id: LocalCrystalPlastFp2D.h,v 1.6 2004-09-10 22:39:43 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_FP_2D_H_
#define _LOCAL_CRYSTAL_PLAST_FP_2D_H_

#include "LocalCrystalPlastFp.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class LocalCrystalPlastFp2D : public LocalCrystalPlastFp
{
 public:
  // constructor
  LocalCrystalPlastFp2D(ifstreamT& in, const FSMatSupportT& support);

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

 protected:
 
  // crystal Cauchy stress in 2D
  dSymMatrixT f2Dsavg_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dcavg_ijkl;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_FP_2D_H_ */
