/* $Id: BCJHypoIsoDamageKE2D.h,v 1.5 2004-07-15 08:29:14 paklein Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_
#define _BCJ_HYPO_ISO_DAMAGE_KE_2D_H_

#include "BCJHypoIsoDamageKE3D.h"

#include <iostream.h>
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

 protected:

  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

} // namespace Tahoe 
#endif /* _BCJ_HYPO_ISO_DAMAGE_KE_2D_ */
