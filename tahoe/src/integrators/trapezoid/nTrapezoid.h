/* $Id: nTrapezoid.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/03/1999)                                          */

#ifndef _N_TRAPEZOID_H_
#define _N_TRAPEZOID_H_

/* base class */
#include "ControllerT.h"
#include "nDtControllerT.h"

/* forward declarations */
class dArrayT;
class KBC_CardT;

class nTrapezoid: public virtual ControllerT, public nDtControllerT
{
public:

	/* constructor */
	nTrapezoid(void);

	/* consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC);
	
	/* predictors - map ALL */
	virtual void Predictor(void);

	/* correctors - map ACTIVE */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update);
	virtual void MappedCorrector(const iArrayT& map, const iArray2DT& flags,
		const dArray2DT& update);

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
