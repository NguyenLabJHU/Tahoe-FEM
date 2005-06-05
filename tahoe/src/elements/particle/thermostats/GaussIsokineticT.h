/* $Id: GaussIsokineticT.h,v 1.4.20.1 2005-06-05 06:23:55 paklein Exp $ */
#ifndef _GAUSS_ISOKINETIC_T_H_
#define _GAUSS_ISOKINETIC_T_H_

/* base class */
#include "ThermostatBaseT.h"

namespace Tahoe {

/** base class for thermostatting and damping */
class GaussIsokineticT: public ThermostatBaseT
{
public:

	/** constructor */
//	GaussIsokineticT(ifstreamT& in, const int& nsd, const double& dt);
	GaussIsokineticT(const BasicSupportT& support);

	/** destructor */
	virtual ~GaussIsokineticT(void) {};
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
};

} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
