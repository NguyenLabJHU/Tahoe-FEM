/* $Id: LocalCrystalPlast2D.h,v 1.5 2003-01-29 07:35:04 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_2D_H_
#define _LOCAL_CRYSTAL_PLAST_2D_H_

#include "LocalCrystalPlast.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class LocalCrystalPlast2D : public LocalCrystalPlast, public Material2DT
{
 public:
  // constructor
  LocalCrystalPlast2D(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~LocalCrystalPlast2D();

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
#endif /* _LOCAL_CRYSTAL_PLAST_2D_H_ */
