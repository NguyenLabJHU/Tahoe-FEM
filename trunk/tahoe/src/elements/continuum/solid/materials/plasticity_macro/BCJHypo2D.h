/* $Id: BCJHypo2D.h,v 1.9 2011-12-01 21:11:38 bcyansfn Exp $ */
#ifndef _BCJ_HYPO_2D_H_
#define _BCJ_HYPO_2D_H_

#include "BCJHypo3D.h"

#include <iostream>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class BCJHypo2D : public BCJHypo3D
{
 public:
  // constructor
  BCJHypo2D(ifstreamT& in, const FSMatSupportT& support);

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

 protected:

  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

} // namespace Tahoe 
#endif /* _BCJ_HYPO_2D_ */
