/* $Id: nTrapezoid.cpp,v 1.2.4.1 2002-04-23 01:24:18 paklein Exp $ */
/* created: paklein (10/03/1999)                                          */

#include "nTrapezoid.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nTrapezoid::nTrapezoid(void): nControllerT(1) { }

/* consistent BC's */
void nTrapezoid::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL)
	{
		cout << "\n nTrapezoid::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU[0])(node, dof);
	double& v = (*fU[1])(node, dof);
	
	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		{
			d = 0.0;
			v = 0.0;
			break;
		}
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			double d_next = KBC.Value();
			v = (d_next - d)/dcorr_v;
			d = d_next;
			break;
		}
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			v  = KBC.Value();
			d += dcorr_v*v;
			break;
		}
		default:
			cout << "\n nTrapezoid::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* predictors - map ALL */
void nTrapezoid::Predictor(void)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL)
	{
		cout << "\n nTrapezoid::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU[0]->Length() != fU[1]->Length()) throw eSizeMismatch;
#endif
	
	/* displacement predictor */
	fU[0]->AddScaled(dpred_v, *fU[1]);
}		

/* correctors - map ACTIVE */
void nTrapezoid::Corrector(const iArray2DT& eqnos, const dArrayT& update,
	int eq_start, int num_eq)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL)
	{
		cout << "\n nTrapezoid::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU[0]->Length()   ||
	   fU[0]->Length() != fU[1]->Length()) throw eSizeMismatch;		
	//NOTE: no check on length of update.
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = fU[0]->Pointer();
	double *pv  = fU[1]->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		/* active dof */
		if (eq > -1 && eq < num_eq)
		{
			double v = update[eq];
			*pd += dcorr_v*v;
			*pv += v;
		}
		pd++;
		pv++;
	}
}

void nTrapezoid::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update, int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL)
	{
		cout << "\n nTrapezoid::MappedCorrector: field arrays not set" << endl;
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
		double* pv = (*fU[1])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag >= eq_start && *pflag <= eq_stop)
			{
				double v = *pupdate;
				*pd += dcorr_v*v;
				*pv += v;
			}
			
			/* next */
			pflag++; pupdate++; pd++; pv++;
		}
	}
}

/* return the field array needed by nControllerT::MappedCorrector. */
const dArray2DT& nTrapezoid::MappedCorrectorField(void) const
{
#if __option (extended_errorcheck)
	if (fU[1] == NULL)
	{
		cout << "\n nTrapezoid::MappedCorrectorField: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	return *fU[1];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nTrapezoid::ExternalNodeCondition(void) const
{
	return KBC_CardT::kVel;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nTrapezoid::nComputeParameters(void)
{
	/* predictor */
	dpred_v = 0.5*fdt;

	/* corrector */
	dcorr_v = 0.5*fdt;
}
