/* $Id: GaussIsokineticT.h,v 1.3 2003-10-30 17:15:21 paklein Exp $ */
#ifndef _GAUSS_ISOKINETIC_T_H_
#define _GAUSS_ISOKINETIC_T_H_

#include "ios_fwd_decl.h"

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** base class for thermostatting and damping */
class GaussIsokineticT: public ThermostatBaseT
{
public:

	/** constructor */
	GaussIsokineticT(ifstreamT& in, const int& nsd, const double& dt);
	GaussIsokineticT(void);

	/** destructor */
	virtual ~GaussIsokineticT(void) {};
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
					
	/** write properties to output */
	virtual void Write(ostream& out) const;
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
protected:

	/** \name properties */
	/*@{*/
	
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
