/* $Id: RampedDampingT.h,v 1.2.12.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _RAMPED_DAMPING_T_H_
#define _RAMPED_DAMPING_T_H_

#include "ios_fwd_decl.h"

#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "RandomNumberT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArray2DT;

/** base class for thermostatting and damping */
class RampedDampingT: public ThermostatBaseT
{
public:

	/** constructor */
	RampedDampingT(ifstreamT& in, const int& nsd, const double& dt);

	/** destructor */
	virtual ~RampedDampingT(void) {};
	
	/** write properties to output */
	virtual void Write(ostream& out) const;
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
	
protected:

	/** \name Damping function variables */
	/*@{*/
	/** values at corners */
	double fLLval, fURval;
	
	/** direction it varies over (only 1D shape function for now)*/
	int nRampedDOF;
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _RAMPED_DAMPING_T_H_ */
