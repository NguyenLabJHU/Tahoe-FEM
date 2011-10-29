/* $Id: LocalCrystalPlast2D.h,v 1.8.28.1 2011-10-29 06:09:09 bcyansfn Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_2D_H_
#define _LOCAL_CRYSTAL_PLAST_2D_H_

#include "LocalCrystalPlast.h"

#include <iostream>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class LocalCrystalPlast2D : public LocalCrystalPlast
{
 public:
  // constructor
  LocalCrystalPlast2D(void);

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

 protected:
 
  // crystal Cauchy stress in 2D
  dSymMatrixT f2Dsavg_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dcavg_ijkl;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_2D_H_ */
