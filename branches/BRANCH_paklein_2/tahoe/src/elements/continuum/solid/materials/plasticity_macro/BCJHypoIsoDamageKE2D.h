/* $Id: BCJHypoIsoDamageKE2D.h,v 1.2.8.1 2002-10-28 06:49:24 paklein Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_
#define _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_

#include "BCJHypoIsoDamageKE3D.h"
#include "Material2DT.h"

#include <iostream.h>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class ElasticT;

class BCJHypoIsoDamageKE2D : public BCJHypoIsoDamageKE3D, public Material2DT
{
 public:
  // constructor
  BCJHypoIsoDamageKE2D(ifstreamT& in, const FDMatSupportT& support);

  // destructor
  ~BCJHypoIsoDamageKE2D();

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:

  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

} // namespace Tahoe 
#endif /* _BCJ_HYPO_ISO_DAMAGE_KE_2D_ */