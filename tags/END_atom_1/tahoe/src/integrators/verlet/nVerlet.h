#ifndef _N_VERLET_H_
#define _N_VERLET_H_

/* base class */
#include "Verlet.h"
#include "nControllerT.h"

namespace Tahoe {

/** Node controller for an explicit 4th order accurate, Verlet time integration
 * algorithm. */
class nVerlet: public virtual Verlet, public nControllerT
{
public:

	/** constructor */
	nVerlet(void);

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