/* $Id: nLinearStaticIntegrator.cpp,v 1.3 2002-07-02 19:55:09 cjkimme Exp $ */
/* created: paklein (10/14/1996) */

#include "nLinearStaticIntegrator.h"
#include "BasicFieldT.h"
#include "dArrayT.h"

/* constructor */

using namespace Tahoe;

nLinearStaticIntegrator::nLinearStaticIntegrator(void) { };

/* predictor. Maps ALL degrees of freedom forward. */
void nLinearStaticIntegrator::Predictor(BasicFieldT& field)
{
	/* clear all displacements */
	field[0] = 0.0;
}

/* correctors - map ACTIVE */
void nLinearStaticIntegrator::Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq)
{
	const iArray2DT& eqnos = field.Equations();

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = field[0].Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
	
		/* active dof */
		if (eq > -1 && eq < num_eq) *pd = update[eq];

		/* next */
		pd++;
	}
}
