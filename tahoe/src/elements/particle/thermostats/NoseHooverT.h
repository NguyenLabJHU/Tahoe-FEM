/* $Id: NoseHooverT.h,v 1.5.30.1 2005-07-25 02:37:14 paklein Exp $ */
#ifndef _NOSE_HOOVER_T_H_
#define _NOSE_HOOVER_T_H_

/* base class */
#include "ThermostatBaseT.h"

namespace Tahoe {

/** base class for thermostatting and damping */
class NoseHooverT: public ThermostatBaseT
{
public:

	/** constructor */
//	NoseHooverT(ifstreamT& in, const int& nsd, const double& dt);
	NoseHooverT(const BasicSupportT& support);

	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
	
	/** \name restart information */
	/*@{*/
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** \name properties */
	/*@{*/
	double fBetaOrig;
	double fEta;
	double fEtaDot;
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
