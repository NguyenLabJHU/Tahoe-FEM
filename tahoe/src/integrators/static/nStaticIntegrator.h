/* $Id: nStaticIntegrator.h,v 1.3.2.1 2002-04-23 01:24:17 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _N_STATIC_CONTROLLER_H_
#define _N_STATIC_CONTROLLER_H_

/* base classes */
#include "StaticT.h"
#include "nIntegratorT.h"

/** nodal integrator for quasistatic systems */
class nStaticIntegrator: public virtual StaticT, public nIntegratorT
{
public:

	/** constructor */
	nStaticIntegrator(void);

	/** consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(void);

	/** corrector - map ACTIVE. See nIntegratorT::Corrector for more
	 * documentation */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update,
		int eq_start, int num_eq);

	/** corrector with node number map - map ACTIVE. See 
	 * nIntegratorT::MappedCorrector for more documentation */
	virtual void MappedCorrector(const iArrayT& map, const iArray2DT& eqnos,
		const dArray2DT& update, int eq_start, int eq_stop);

	/** return the field array needed by nIntegratorT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(void) const;
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void nComputeParameters(void);	
};

#endif /* _N_STATIC_CONTROLLER_H_ */
