/* $Id: nLinearStaticIntegrator.cpp,v 1.1.4.1 2002-04-24 01:29:22 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "nLinearStaticIntegrator.h"
#include "BasicFieldT.h"

/* constructor */
nLinearStaticIntegrator::nLinearStaticIntegrator(void) { };

/* predictor. Maps ALL degrees of freedom forward. */
void nLinearStaticIntegrator::Predictor(BasicFieldT& field)
{
	/* clear all displacements */
	field[0] = 0.0;
}
