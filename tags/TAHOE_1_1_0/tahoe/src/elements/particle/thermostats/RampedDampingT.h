/* $Id: RampedDampingT.h,v 1.1 2003-04-18 19:01:56 cjkimme Exp $ */
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
	RampedDampingT(ifstreamT& in, int nsd, double dt);

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

	/** \name properties */
	/*@{*/
	double fBeta;	
	double fTemperature;
	double fTimeStep;
	/*@}*/
	
	/** Nodes that are thermostatted */
	iArrayT fNodes;
	
	/** True if stochastic force is returned */
//	bool QLangevin;
	
	/** Number of spatial dimensions */
	int fSD;
	
//	RandomNumberT* fRandom;
	bool qNodesInRegion;
	
	/** \name Region paramters */
	/*@{*/
	/** Bounding box */
	dArrayT fxmin, fxmax;
	/** Re-check for particles in region every nIncs timesteps */
	int nIncs; 
	/*@}*/
};


} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
