/* $Id: nNLHHTalpha.h,v 1.2 2001-08-27 17:12:11 paklein Exp $ */
/* created: paklein (10/17/1996) */

#ifndef _N_NL_HHT_A_H_
#define _N_NL_HHT_A_H_

/* base classes */
#include "HHTalpha.h"
#include "nControllerT.h"

/** HHT alpha integration for linear systems */
class nNLHHTalpha: public virtual HHTalpha, public nControllerT
{
public:

	/* constructor */
	nNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder = kHHTalphaAuto_O2);

	/* consistent BC's - updates predictors and acceleration only */
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
	double	dpred_a;
	double	vpred_a;
	
	/* corrector/consistent BC */
	double	dcorr_a;
	double	vcorr_a;
	
};

#endif /* _N_NL_HHT_A_H_ */
