/* $Id: nLinearStaticIntegrator.h,v 1.4 2002-07-05 22:27:55 paklein Exp $ */
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

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(BasicFieldT& field);	  	

	/** corrector - map ACTIVE. See nIntegratorT::Corrector for more
	 * documentation */
	virtual void Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq);
};

} // namespace Tahoe 
#endif /* _N_LINEAR_STATIC_CONTROLLER_H_ */