/* $Id: nTrapezoid.h,v 1.8 2003-01-29 07:35:18 paklein Exp $ */
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

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(BasicFieldT& field);

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update);

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
