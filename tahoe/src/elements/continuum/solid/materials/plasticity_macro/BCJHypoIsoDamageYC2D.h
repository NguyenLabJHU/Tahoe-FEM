/* $Id: BCJHypoIsoDamageYC2D.h,v 1.6 2004-09-10 22:39:48 paklein Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_YC_2D_H_
#define _BCJ_HYPO_ISO_DAMAGE_YC_2D_H_

#include "BCJHypoIsoDamageYC3D.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class BCJHypoIsoDamageYC2D : public BCJHypoIsoDamageYC3D
{
 public:
  // constructor
  BCJHypoIsoDamageYC2D(ifstreamT& in, const FSMatSupportT& support);

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
#endif /* _BCJ_HYPO_ISO_DAMAGE_YC_2D_ */
