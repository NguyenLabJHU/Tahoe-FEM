/* $Id: nTrapezoid.h,v 1.8.60.2 2004-11-15 04:15:00 d-farrell2 Exp $ */
/* created: paklein (10/03/1999) */
#ifndef _N_TRAPEZOID_H_
#define _N_TRAPEZOID_H_

/* base class */
#include "Trapezoid.h"
#include "nIntegratorT.h"

namespace Tahoe {

/** trapezoidal integration for first order systems */
class nTrapezoid: public virtual Trapezoid, public nIntegratorT
{
public:

	/* constructor */
	nTrapezoid(void);

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
	double	dpred_v;
	
	/* corrector */
	double	dcorr_v;
};

} // namespace Tahoe 
#endif /* _N_TRAPEZOID_H_ */
