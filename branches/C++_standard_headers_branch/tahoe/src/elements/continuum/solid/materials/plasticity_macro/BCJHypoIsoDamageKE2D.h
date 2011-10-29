/* $Id: BCJHypoIsoDamageKE2D.h,v 1.6.36.1 2011-10-29 06:09:09 bcyansfn Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_
#define _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_

#include "BCJHypoIsoDamageKE3D.h"

#include <iostream>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class BCJHypoIsoDamageKE2D : public BCJHypoIsoDamageKE3D
{
 public:
  // constructor
  BCJHypoIsoDamageKE2D(ifstreamT& in, const FSMatSupportT& support);

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
#endif /* _BCJ_HYPO_ISO_DAMAGE_KE_2D_ */
