/* $Id: nStaticIntegrator.cpp,v 1.2.4.1 2002-04-23 01:24:17 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "nStaticIntegrator.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* constructor */
nStaticIntegrator::nStaticIntegrator(void):nIntegratorT(0) { }

/* consistent BC's */
void nStaticIntegrator::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nStaticIntegrator::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destination */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU[0])(node, dof);
	
	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		{
			d = 0.0;
			break;
		}
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			d = KBC.Value();
			break;
		}
		default:
			cout << "\n nTrapezoid::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* predictor. Maps ALL degrees of freedom forward. */
void nStaticIntegrator::Predictor(void)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nStaticIntegrator::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	//nothing to do
}

/* correctors - map ACTIVE */
void nStaticIntegrator::Corrector(const iArray2DT& eqnos, const dArrayT& update,
	int eq_start, int num_eq)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nStaticIntegrator::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU[0]->Length()) throw eSizeMismatch;		
	//NOTE: no check on length of update.
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = fU[0]->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
	
		/* active dof */
		if (eq > -1 && eq < num_eq) *pd += update[eq];

		/* next */
		pd++;
	}
}

void nStaticIntegrator::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update, int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nStaticIntegrator::MappedCorrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* checks */
	if (flags.MajorDim() != map.Length() ||
	    flags.MajorDim() != update.MajorDim() ||
	    flags.MinorDim() != update.MinorDim() ||
	    flags.MinorDim() != fU[0]->MinorDim()) throw eSizeMismatch;

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		int* pflags = flags(i);

		double* pd = (*fU[0])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag >= eq_start && *pflag <= eq_stop) *pd = *pupdate;
			
			/* next */
			pflag++; pupdate++; pd++;
		}
	}
}

/* return the field array needed by nIntegrator::MappedCorrector. */
const dArray2DT& nStaticIntegrator::MappedCorrectorField(void) const
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL)
	{
		cout << "\n nExplicitCD::MappedCorrectorField: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	return *fU[0];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nStaticIntegrator::ExternalNodeCondition(void) const
{
	return KBC_CardT::kDsp;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate time stepping constants */
void nStaticIntegrator::nComputeParameters(void) { }
