/* $Id: LangevinT.h,v 1.4.22.1 2004-05-25 16:36:43 paklein Exp $ */
#ifndef _LANGEVIN_T_H_
#define _LANGEVIN_T_H_

/* base class */
#include "ThermostatBaseT.h"

namespace Tahoe {

/* forward declarations */
class RandomNumberT;

/** insert witty comment here */
class LangevinT: public ThermostatBaseT
{
public:
	
	/** constructor */
//	LangevinT(ifstreamT& in, const int& nsd, const double& dt);
	LangevinT(const BasicSupportT& support);

	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
	
	void SetRandNumGenerator(RandomNumberT* frand);
	
protected:

	/** \name properties */
	/*@{*/
//	double fAmp;
	/*@}*/
	
	RandomNumberT* fRandom;
};

inline void LangevinT::SetRandNumGenerator(RandomNumberT* frand)
{
	fRandom = frand;
}


} /* namespace Tahoe */

#endif /* _LANGEVIN_T_H_ */
