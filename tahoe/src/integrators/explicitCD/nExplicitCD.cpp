/* $Id: nExplicitCD.cpp,v 1.2.4.1 2002-04-23 01:24:16 paklein Exp $ */
/* created: paklein (03/23/1997) */

#include "nExplicitCD.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nExplicitCD::nExplicitCD(void):nControllerT(2) { }

/* consistent BC's - updates predictors and acceleration only */
void nExplicitCD::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nExplicitCD::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU[0])(node, dof);
	double& v = (*fU[1])(node, dof);
	double& a = (*fU[2])(node, dof);

	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		{
			d = 0.0;
			v = 0.0; //correct?	
			a = 0.0; //correct?	
			break;
		}
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			d = KBC.Value();
			a = 0.0; //correct?	
			break;
		}
		
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			double v_next = KBC.Value();
	
			a = (v_next - v)/vcorr_a;
			v = v_next;
			break;
		}
		
		case KBC_CardT::kAcc: /* prescribed acceleration */
		{
			a  = KBC.Value();
			v += vcorr_a*a;
			break;
		}

		default:
		
			cout << "\n nExplicitCD::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* predictors - map ALL */
void nExplicitCD::Predictor(void)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nExplicitCD::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU[0]->Length() != fU[1]->Length() ||
	    fU[1]->Length() != fU[2]->Length()) throw eGeneralFail;
#endif

	/* displacement predictor */
	fU[0]->AddCombination(dpred_v, *fU[1], dpred_a, *fU[2]);	

	/* velocity predictor */
	fU[1]->AddScaled(vpred_a, *fU[2]);
}		

/* correctors - map ACTIVE */
void nExplicitCD::Corrector(const iArray2DT& eqnos, const dArrayT& update,
	int eq_start, int num_eq)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nExplicitCD::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU[0]->Length() ||
	   fU[0]->Length() != fU[1]->Length()  ||
	   fU[1]->Length() != fU[2]->Length()) throw eGeneralFail;		
	//NOTE: no check on length of update.
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pv  = fU[1]->Pointer();
	double *pa  = fU[2]->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		
		/* active dof */
		if (eq > -1 && eq < num_eq)
		{
			double a = update[eq];
			*pv += vcorr_a*a;
			*pa = a;
		}
		pv++;
		pa++;
	}
}

void nExplicitCD::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update, int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nExplicitCD::MappedCorrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* checks */
	if (flags.MajorDim() != map.Length() ||
	    flags.MajorDim() != update.MajorDim() ||
	    flags.MinorDim() != update.MinorDim() ||
	    flags.MinorDim() != fU[1]->MinorDim()) throw eSizeMismatch;

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		int* pflags = flags(i);

		double* pv = (*fU[1])(row);
		double* pa = (*fU[2])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag >= eq_start && *pflag <= eq_stop)
			{
				double a = *pupdate;
				*pv += vcorr_a*a;
				*pa = a;
			}
			
			/* next */
			pflag++; pupdate++; pv++; pa++;
		}
	}
}

/* return the field array needed by nControllerT::MappedCorrector. */
const dArray2DT& nExplicitCD::MappedCorrectorField(void) const
{
#if __option (extended_errorcheck)
	if (fU[2] == NULL)
	{
		cout << "\n nExplicitCD::MappedCorrectorField: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	return *fU[2];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nExplicitCD::ExternalNodeCondition(void) const
{
	return KBC_CardT::kAcc;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nExplicitCD::nComputeParameters(void)
{
	/* predictor */
	dpred_v		= fdt;
	dpred_a		= 0.5*fdt*fdt;
	vpred_a		= 0.5*fdt;
	
	/* corrector */
	vcorr_a		= 0.5*fdt;
}
