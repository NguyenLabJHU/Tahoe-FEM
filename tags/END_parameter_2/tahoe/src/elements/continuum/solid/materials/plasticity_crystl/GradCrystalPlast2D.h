/* $Id: GradCrystalPlast2D.h,v 1.5.30.1 2004-03-03 16:15:04 paklein Exp $ */
#ifndef _GRAD_CRYSTAL_PLAST_2D_H_
#define _GRAD_CRYSTAL_PLAST_2D_H_

#include "GradCrystalPlast.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class GradCrystalPlast2D: public GradCrystalPlast
{
 public:
  // constructor
  GradCrystalPlast2D(ifstreamT& in, const FSMatSupportT& support);

  // crystal Cauchy stress
  virtual const dSymMatrixT& s_ij();

  // crystal modulus 
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void PrintName(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

 protected: 

  // crystal Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dc_ijkl;

};

} // namespace Tahoe 
#endif /* _GRAD_CRYSTAL_PLAST_2D_H_ */
