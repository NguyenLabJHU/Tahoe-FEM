/* $Id: nTrapezoid.h,v 1.2 2001-08-27 17:12:17 paklein Exp $ */
/* created: paklein (10/03/1999) */

#ifndef _N_TRAPEZOID_H_
#define _N_TRAPEZOID_H_

/* base class */
#include "ControllerT.h"
#include "nControllerT.h"

/** trapezoidal integration for first order systems */
class nTrapezoid: public virtual ControllerT, public nControllerT
{
public:

	/* constructor */
	nTrapezoid(void);

	/* consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC);
	
	/* predictors - map ALL */
	virtual void Predictor(void);

	/** corrector - map ACTIVE. See nControllerT::Corrector for more
	 * documentation */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update,
		int eq_start, int eq_stop);

	/** corrector with node number map - map ACTIVE. See 
	 * nControllerT::MappedCorrector for more documentation */
	virtual void MappedCorrector(const iArrayT& map, const iArray2DT& eqnos,
		const dArray2DT& update, int eq_start, int eq_stop);

	/** return the field array needed by nControllerT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(void) const;

	/* pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void);

private:

	/* predictor */
	double	dpred_v;
	
	/* corrector */
	double	dcorr_v;
};

#endif /* _N_TRAPEZOID_H_ */
