/* $Id: nLinearHHTalpha.h,v 1.6.10.1 2002-12-18 09:38:58 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _N_LINEARHHT_A_H_
#define _N_LINEARHHT_A_H_

/* base class */
#include "HHTalpha.h"
#include "nControllerT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** HHT alpha integration for linear systems */
class nLinearHHTalpha: virtual public HHTalpha, public nControllerT
{
public:

	/** constructor */
	nLinearHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(BasicFieldT& field);

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update);

	/** corrector - map ACTIVE. See nControllerT::Corrector for more
	 * documentation */
	virtual void Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq);

	/** corrector with node number map - map ACTIVE. See 
	 * nControllerT::MappedCorrector for more documentation */
	virtual void MappedCorrector(BasicFieldT& field, const iArrayT& map,
		const iArray2DT& flags, const dArray2DT& update);

	/** return the field array needed by nControllerT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(BasicFieldT& field) const;

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

} // namespace Tahoe 
#endif /* _N_LINEARHHT_A_H_ */
