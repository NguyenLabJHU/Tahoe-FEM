/* $Id: nVerlet.h,v 1.7 2004-12-26 21:09:26 d-farrell2 Exp $ */
#ifndef _N_VERLET_H_
#define _N_VERLET_H_

/* base class */
#include "Verlet.h"
#include "nIntegratorT.h"

namespace Tahoe {

/** Node controller for an explicit 4th order accurate, Verlet time integration
 * algorithm. */
class nVerlet: public virtual Verlet, public nIntegratorT
{
public:

	/** constructor */
	nVerlet(void);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward, Unless specified otherwise*/
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);

	/** corrector. Maps ALL degrees of freedom forward, Unless specified otherwise  */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart = 0, int fieldend = -1, int dummy = 0);

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
	
	/** recalculate time stepping constants */
	virtual void nComputeParameters(void);
	
private:

	/** \name predictor constants*/
	/*@{*/
	double dpred_v;
	double dpred_a;
	double vpred_a;
	/*@}*/
	
	/** corrector constant, also used for consistent BC*/  	
	double vcorr_a;
	  	  	
};

} // namespace Tahoe

#endif /* _N_VERLET_H_ */
