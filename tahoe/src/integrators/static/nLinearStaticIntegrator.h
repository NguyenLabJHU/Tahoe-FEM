/* $Id: nLinearStaticIntegrator.h,v 1.5.60.1 2004-11-08 02:16:02 d-farrell2 Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _N_LINEAR_STATIC_CONTROLLER_H_
#define _N_LINEAR_STATIC_CONTROLLER_H_

/* base class */
#include "nStaticIntegrator.h"

namespace Tahoe {

/** nodal integrator for linear quasistatic systems. This
 * integrator differs from nStaticIntegratorT only in that
 * the predictor for linear systems sets the displacement
 * of every node to 0.0. */
class nLinearStaticIntegrator: public nStaticIntegrator
{
public:

	/** constructor */
	nLinearStaticIntegrator(void);

	/** predictor. Maps ALL degrees of freedom forward, Unless specified otherwise */
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);	  	

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update);

	/** corrector - map ACTIVE. See nIntegratorT::Corrector for more
	 * documentation */
	virtual void Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq);
};

} // namespace Tahoe 
#endif /* _N_LINEAR_STATIC_CONTROLLER_H_ */
