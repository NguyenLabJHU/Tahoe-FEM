/* $Id: LocalCrystalPlastFp2D.h,v 1.3 2002-11-14 17:06:32 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_FP_2D_H_
#define _LOCAL_CRYSTAL_PLAST_FP_2D_H_

#include "LocalCrystalPlastFp.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class ElasticT;

class LocalCrystalPlastFp2D : public LocalCrystalPlastFp, public Material2DT
{
 public:
  // constructor
  LocalCrystalPlastFp2D(ifstreamT& in, const FDMatSupportT& support);

  // destructor
  ~LocalCrystalPlastFp2D();

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:
 
  // crystal Cauchy stress in 2D
  dSymMatrixT f2Dsavg_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dcavg_ijkl;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_FP_2D_H_ */
