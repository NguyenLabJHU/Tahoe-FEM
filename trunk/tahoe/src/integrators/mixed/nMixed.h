/* $Id: */
/* created: a-kopacz (08/08/2006) */

#ifndef _N_MIXED_H_
#define _N_MIXED_H_

/* base class */
#include "Mixed.h"
#include "nIntegratorT.h"
#include "dArrayT.h"

namespace Tahoe {

/** mixed integration for first order systems */
class nMixed: public virtual Mixed, public nIntegratorT
{
public:

	/* constructor */
	nMixed(void);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward Unless specified otherwise */
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);

	/** corrector. Maps ALL degrees of freedom forward Unless specified otherwise*/
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart = 0, int fieldend = -1, int dummy = 0);

	/** corrector - map ACTIVE. See nIntegratorT::Corrector for more
	 * documentation */
	virtual void Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq);

	/** corrector with node number map - map ACTIVE. See 
	 * nIntegratorT::MappedCorrector for more documentation */
	virtual void MappedCorrector(BasicFieldT& field, const iArrayT& map, 
		const iArray2DT& flags, const dArray2DT& update);

	/** return the field array needed by nIntegratorT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(BasicFieldT& field) const;

protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void);

private:

	/* predictor */
	dArrayT	dpred_v_;
	
	/* corrector */
	dArrayT	dcorr_v_;
};

} // namespace Tahoe 
#endif /* _N_MIXED_H_ */
