/* $Id: nLinearHHTalpha.h,v 1.3.2.1 2002-04-23 01:24:15 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _N_LINEARHHT_A_H_
#define _N_LINEARHHT_A_H_

/* base class */
#include "HHTalpha.h"
#include "nControllerT.h"

/* direct members */
#include "dArray2DT.h"

/** HHT alpha integration for linear systems */
class nLinearHHTalpha: virtual public HHTalpha, public nControllerT
{
public:

	/** constructor */
	nLinearHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder);

	/* consistent BC's - updates predictors and acceleration only */
	virtual void ConsistentKBC(const KBC_CardT& KBC);
	
	/* predictors - map ALL */
	virtual void Predictor(void);

	/** corrector - map ACTIVE. See nControllerT::Corrector for more
	 * documentation */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update,
		int eq_start, int num_eq);

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

	/* values from t_n (needed for HHT-alpha) */
	dArray2DT	dn;
	dArray2DT	vn;

	/* predictor */
	double	dpred_v;
	double	dpred_a;
	double	vpred_a;
	
	/* corrector */
	double	dcorr_d;
	double	dcorr_dpred;
	double	dcorr_a;
	
	double	vcorr_v;
	double	vcorr_vpred;
	double	vcorr_a;
	  	
	/* consistent BC */
	double	dalpha_a;
	double	valpha_a;
		
};

#endif /* _N_LINEARHHT_A_H_ */
