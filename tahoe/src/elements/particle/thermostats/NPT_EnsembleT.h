/* $Id: NPT_EnsembleT.h,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _NPT_ENSEMBLE_T_H_
#define _NPT_ENSEMBLE_T_H_

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

/** base class for thermostatting and damping */
class NPT_EnsembleT: public ThermostatBaseT
{
public:

	/** constructor */
	NPT_EnsembleT(ifstreamT& in, const ElementSupportT& support,
		const dArray2DT& dynStress);

	/** destructor */
	virtual ~NPT_EnsembleT(void) {};
	
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
	const ElementSupportT& fElementSupport;
	const dArray2DT& fDynStress;
	/** Needed to control periodic boundary lengths as volume changes */
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _NPT_ENSEMBLE_T_H_ */
