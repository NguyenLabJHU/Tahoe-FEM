/* $Id: NoseHooverT.h,v 1.3.12.1 2003-09-18 21:03:38 cjkimme Exp $ */
#ifndef _NOSE_HOOVER_T_H_
#define _NOSE_HOOVER_T_H_

#include "ios_fwd_decl.h"

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Feedback method for thermostatting. See PRA _31_ 1695 for details.  */
class NoseHooverT: public ThermostatBaseT
{
public:

	/** constructor */
	NoseHooverT(ifstreamT& in, const int& nsd, const double& dt);

	/** destructor */
	virtual ~NoseHooverT(void) {};
	
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
	/*@}*/
	
};

} /* namespace Tahoe */

#endif /* _NOSE_HOOVER_T_H_ */
