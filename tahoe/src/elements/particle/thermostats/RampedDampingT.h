/* $Id: RampedDampingT.h,v 1.4.20.1 2005-06-05 06:23:55 paklein Exp $ */
#ifndef _RAMPED_DAMPING_T_H_
#define _RAMPED_DAMPING_T_H_

/* base class */
#include "ThermostatBaseT.h"

namespace Tahoe {

/** base class for thermostatting and damping */
class RampedDampingT: public ThermostatBaseT
{
public:

	/** constructor */
//	RampedDampingT(ifstreamT& in, const int& nsd, const double& dt);
	RampedDampingT(const BasicSupportT& support);
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
	
protected:

	/** \name properties */
	/*@{*/
	double fBeta;
	/*@}*/

	bool qNodesInRegion;	
};

} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
