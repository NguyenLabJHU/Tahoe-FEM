/* $Id: RampedDampingT.h,v 1.4 2004-07-15 08:29:54 paklein Exp $ */
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
			dArray2DT& forces, AutoArrayT<int>& types,
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
