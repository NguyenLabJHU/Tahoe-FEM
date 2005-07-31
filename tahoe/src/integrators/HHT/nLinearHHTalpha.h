/* $Id: nLinearHHTalpha.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#ifndef _N_LINEARHHT_A_H_
#define _N_LINEARHHT_A_H_

/* base class */
#include "HHTalpha.h"
#include "nDDtControllerT.h"

/* direct members */
#include "dArray2DT.h"

/* forward declarations */
class dArrayT;
class KBC_CardT;

class nLinearHHTalpha: public virtual HHTalpha, public nDDtControllerT
{
public:

	/* constructor */
	nLinearHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder = kHHTalphaAuto_O2);

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
