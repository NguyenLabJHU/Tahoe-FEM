/* $Id: ConstrainedPressureT.h,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _CONSTRAINED_PRESSURE_T_H_
#define _CONSTRAINED_PRESSURE_T_H_

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

/** Constraint method for constant (virial) pressure simulations.
    See Phys. Lett. _98A_ 433 and Comput. Phys. Rep. _1_ 297
    for details. */
class ConstrainedPressureT: public ThermostatBaseT
{
public:

	/** constructor */
	ConstrainedPressureT(ifstreamT& in, const ElementSupportT& support,
		const dArray2DT& dynStress);

	/** destructor */
	virtual ~ConstrainedPressureT(void) {};
	
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
	
	static double VelocityCorrector(void);
	
protected:

	/** \name properties */
	/*@{*/
	const ElementSupportT& fElementSupport;
	const dArray2DT& fDynStress;
	/** Needed to control periodic boundary lengths as volume changes */
	/*@}*/
	
	static double staticChi;
	/** velocity corrector needed by KBCController */
	
	double fV; 
	/** current simulation volume */
	double fVdot;
	/** rate of change of simulation volume */
	
};

} /* namespace Tahoe */

#endif /* _CONSTRAINED_PRESSURE_T_H_ */
