/* $Id: ConfigurationalT.h,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _CONFIGURATIONAL_T_H_
#define _CONFIGURATIONAL_T_H_

#include "ios_fwd_decl.h"

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ElementSupportT;

/** base class for thermostatting and damping */
class ConfigurationalT: public ThermostatBaseT
{
public:

	/** constructor */
	ConfigurationalT(ifstreamT& in, const ElementSupportT& support, 
		const dArrayT& delDotF);

	/** destructor */
	virtual ~ConfigurationalT(void) {};
	
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
	double fBetaOrig;
	double fEta;
	double fEtaDot;
	const dArrayT& fDelDotF;
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _CONFIGURATIONAL_T_H_ */
