/* $Id: nNLHHTalpha.h,v 1.10.4.1 2004-11-08 02:15:57 d-farrell2 Exp $ */
/* created: paklein (10/17/1996) */
#ifndef _N_NL_HHT_A_H_
#define _N_NL_HHT_A_H_

/* base classes */
#include "HHTalpha.h"
#include "nIntegratorT.h"

namespace Tahoe {

/** HHT alpha integration for linear systems */
class nNLHHTalpha: public virtual HHTalpha, public nIntegratorT
{
public:

	/** constructor */
	nNLHHTalpha(double alpha);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward Unless specified otherwise */
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);

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
	double	dpred_a;
	double	vpred_a;
	
	/* corrector/consistent BC */
	double	dcorr_a;
	double	vcorr_a;
};

} // namespace Tahoe 
#endif /* _N_NL_HHT_A_H_ */
