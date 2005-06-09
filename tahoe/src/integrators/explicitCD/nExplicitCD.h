/* $Id: nExplicitCD.h,v 1.9.12.2 2005-06-09 02:43:56 d-farrell2 Exp $ */
/* created: paklein (03/23/1997) */

#ifndef _N_EXP_CD_H_
#define _N_EXP_CD_H_

/* base class */
#include "ExplicitCD.h"
#include "nIntegratorT.h"

namespace Tahoe {

/** Node controller for an explicit 2nd order accurate, central 
 * difference time-stepping algorithm. */
class nExplicitCD: public virtual ExplicitCD, public nIntegratorT
{
public:

	/** constructor */
	nExplicitCD(void);

	/** prescribe the field and derivatives consistent BC's 
	 * \param field field to which boundary conditions will be applied
	 * \param KBC boundary condition specification
	 * \param nodes pointer to a list of nodes if KBC card has KBC_CardT::ModeT set to KBC_CardT::kSet
	 *        or NULL if it is set to KBC_CardT::kNode. */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC, const iArrayT* nodes = NULL);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward Unless specified otherwise */
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart = 0, int fieldend = -1);

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
#endif /* _N_EXP_CD_H_ */
