/* $Id: nLinearStaticIntegrator.h,v 1.1.4.1 2002-04-24 01:29:22 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _N_LINEAR_STATIC_CONTROLLER_H_
#define _N_LINEAR_STATIC_CONTROLLER_H_

/* base class */
#include "nStaticIntegrator.h"

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
};

#endif /* _N_LINEAR_STATIC_CONTROLLER_H_ */
