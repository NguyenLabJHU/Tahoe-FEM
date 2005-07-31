/* $Id: nExplicitCD.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Node controller for an explicit 2nd order                              */
/* accurate, central difference time-stepping algorithm.                  */

#ifndef _N_EXP_CD_H_
#define _N_EXP_CD_H_

/* base class */
#include "ControllerT.h"
#include "nDDtControllerT.h"

class nExplicitCD: public virtual ControllerT, public nDDtControllerT
{
public:

	/* constructor */
	nExplicitCD(void);

	/* consistent BC's - updates predictors and acceleration only */
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
	double	dpred_a;
	double	vpred_a;
	
	/* corrector */  	
	double	vcorr_a; //also used for consistent BC
	  	  	
};

#endif /* _N_EXP_CD_H_ */
