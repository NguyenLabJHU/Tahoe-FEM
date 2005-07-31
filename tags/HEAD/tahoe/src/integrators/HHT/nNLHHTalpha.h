/* $Id: nNLHHTalpha.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/17/1996)                                          */

#ifndef _N_NL_HHT_A_H_
#define _N_NL_HHT_A_H_

/* base classes */
#include "HHTalpha.h"
#include "nDDtControllerT.h"

/* forward declarations */
class dArrayT;
class dArray2DT;
class KBC_CardT;

class nNLHHTalpha: public virtual HHTalpha, public nDDtControllerT
{
public:

	/* constructor */
	nNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder = kHHTalphaAuto_O2);

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
	
	/* corrector/consistent BC */
	double	dcorr_a;
	double	vcorr_a;
	
};

#endif /* _N_NL_HHT_A_H_ */
