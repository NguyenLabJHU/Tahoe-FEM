/* $Id: LangevinT.h,v 1.4.22.2 2004-05-26 03:56:15 paklein Exp $ */
#ifndef _LANGEVIN_T_H_
#define _LANGEVIN_T_H_

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "RandomNumberT.h"

namespace Tahoe {

/** insert witty comment here */
class LangevinT: public ThermostatBaseT
{
public:
	
	/** constructor */
	LangevinT(const BasicSupportT& support);

	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/* random number generator */
	RandomNumberT fRandom;
};

} /* namespace Tahoe */

#endif /* _LANGEVIN_T_H_ */
