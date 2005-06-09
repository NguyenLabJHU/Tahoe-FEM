/* $Id: nGear6.h,v 1.9.12.2 2005-06-09 02:43:57 d-farrell2 Exp $ */
#ifndef _N_GEAR_06_H_
#define _N_GEAR_06_H_

/* base class */
#include "Gear6.h"
#include "nIntegratorT.h"
#include "dArray2DT.h"

namespace Tahoe {

/** Node controller for an explicit 6th order accurate Gear time integration
 * algorithm. The BasicFieldT used with this integrator must have
 * BasicFieldT::Order of 6. */
class nGear6: public virtual Gear6, public nIntegratorT
{
public:

	/** constructor */
	nGear6(void);

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

private:

	/** \name Gear constants */
	/*@{*/
	double F02; 
	double F12;
	double F22;
	double F32;
	double F42;
	double F52;
	/*@}*/
	
	/** \name Taylor expansion factors */
	/*@{*/
	double fdt2; /**< \f$ \frac{\Delta t^2}{2!} \f$ */
	double fdt3; /**< \f$ \frac{\Delta t^3}{3!} \f$ */
	double fdt4; /**< \f$ \frac{\Delta t^4}{4!} \f$ */
	double fdt5; /**< \f$ \frac{\Delta t^5}{5!} \f$ */
	/*@}*/

 protected:

	/** recalculate time stepping constants */
	virtual void nComputeParameters(void);

};

} // namespace Tahoe

#endif /* _N_GEAR_06_H_ */
