/* $Id: AndersenPressureT.h,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _ANDERSEN_PRESSURE_T_H_
#define _ANDERSEN_PRESSURE_T_H_

#include "ios_fwd_decl.h"

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "ElementSupportT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Extended-system method for constant (virial) pressure simulations.
    See JCP _72_ 2384 for details. This method _must_ be used with
    the Gear4 integrator. The parameter Beta is the mass of the
    piston producing the uniform volume changes. */
class AndersenPressureT: public ThermostatBaseT
{
public:

	/** constructor */
	AndersenPressureT(ifstreamT& in, const ElementSupportT& support,
		const double* virial, const double boxLength);

	/** destructor */
	virtual ~AndersenPressureT(void) {};
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
			
	/** close the current timestep. Adjust the system volume */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
					
	/** write properties to output */
	virtual void Write(ostream& out) const;
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
	static double VelocityCorrector(void);
	
protected:

	/** \name properties */
	/*@{*/
	const ElementSupportT& fElementSupport;
	const double* fVirial;
	/** Needed to control periodic boundary lengths as volume changes */
	/*@}*/
	
	dArrayT fV_local;
	/** allocated storage simulation volume and its deriviates */

public:

	static double MdP;
	/** velocity corrector needed by Gear4 integrator */
	
	static dArrayT* fV_field; 
	/** pointer to current simulation volume and its derivatives */
	
};

} /* namespace Tahoe */

#endif /* _CONSTRAINED_PRESSURE_T_H_ */
