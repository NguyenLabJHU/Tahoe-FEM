#ifndef _N_GEAR6_H_
#define _N_GEAR6_H_

/* base class */
#include "Gear6.h"
#include "nControllerT.h"

namespace Tahoe {

/** Node controller for an explicit 6th order accurate Gear time integration
 * algorithm. */
class nGear6: public virtual Gear6, public nControllerT
{
public:

	/** constructor */
	nGear6(void);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(BasicFieldT& field);

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

	/* \name predictor constants*/
	double dpred_v;
	double dpred_a;
	double vpred_a;
	
	/* corrector constant, also used for consistent BC*/  	
	double vcorr_a;

	/* Gear constants */
	double F02; 
	double F12;
	double F32;
	double F42;
	double F52;

	/* higher order derivatives */
	//dArray2DT fD3;
	//dArray2DT fD4; 
	//dArray2DT fD5;
	  	  	
};

} // namespace Tahoe

#endif /* _N_GEAR6_H_ */
