/* $Id: nNLHHTalpha.cpp,v 1.2 2001-08-27 17:12:11 paklein Exp $ */
/* created: paklein (10/17/1996)                                          */

#include "nNLHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nNLHHTalpha::nNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder):
	HHTalpha(in, out, auto2ndorder),
	nControllerT(2)
{

}

/* consistent BC's - updates predictors and acceleration only */
void nNLHHTalpha::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU[0] || !fU[1] || !fU[2])
	{
		cout << "\n nNLHHTalpha::ConsistentKBC: field arrays not set" << endl;
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
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			double temp = KBC.Value();
			
			a  = (temp - d)/dcorr_a;
			d  = temp;
			v += vcorr_a*a;
		
			break;
		}
		
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			double temp = KBC.Value();
			
			a  = (temp - v)/vcorr_a;
			d += dcorr_a*a;
			v  = temp;
	
			break;
		}
		
		case KBC_CardT::kAcc: /* prescribed acceleration */
		{
			a  = KBC.Value();
			d += dcorr_a*a;
			v += vcorr_a*a;

			break;
		}
		
		default:
		
			cout << "\nnNLHHTalpha::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* rredictors - map ALL */
void nNLHHTalpha::Predictor(void)
{
#if __option (extended_errorcheck)
	if (!fU[0] || !fU[1] || !fU[2])
	{
		cout << "\n nNLHHTalpha::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU[0]->Length() != fU[1]->Length() ||
        fU[1]->Length() != fU[2]->Length()) throw eGeneralFail;
#endif

	/* displacement predictor */
	fU[0]->AddCombination(dpred_v, *fU[1], dpred_a, *fU[2]);

	/* velocity predictor */
	fU[1]->AddScaled(vpred_a, *fU[2]);
	
	/* acceleration predictor */
	(*fU[2]) = 0.0;
}		

/* correctors - map ACTIVE */
void nNLHHTalpha::Corrector(const iArray2DT& eqnos, const dArrayT& update,
	int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (!fU[0] || !fU[1] || !fU[2])
	{
		cout << "\n nNLHHTalpha::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU[0]->Length()  ||
       fU[0]->Length() != fU[1]->Length() ||
       fU[1]->Length() != fU[2]->Length()) throw eGeneralFail;		
	/* note: No check on length of update. */
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	
	double *pd  = fU[0]->Pointer();
	double *pv  = fU[1]->Pointer();
	double *pa  = fU[2]->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		
		if (eq > -1 && eq < eq_stop) /* active dof */
		{
			double da = update[eq];
		
			*pd += dcorr_a*da;
			*pv += vcorr_a*da;
			*pa += da;
		}
			
		pd++;
		pv++;
		pa++;		
	}	
}

void nNLHHTalpha::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update, int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (!fU[0] || !fU[1] || !fU[2])
	{
		cout << "\n nNLHHTalpha::MappedCorrector: field arrays not set" << endl;
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
		double* pa = (*fU[2])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag >= eq_start && *pflag <= eq_stop)
			{
				double da = *pupdate;
			
				*pd += dcorr_a*da;
				*pv += vcorr_a*da;
				*pa += da;
			}
			
			/* next */
			pflag++; pupdate++; pd++; pv++; pa++;
		}
	}
}

/* return the field array needed by nControllerT::MappedCorrector. */
const dArray2DT& nNLHHTalpha::MappedCorrectorField(void) const
{
#if __option (extended_errorcheck)
	if (fU[2] == NULL)
	{
		cout << "\n nNLHHTalpha::MappedCorrectorField: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	return *fU[2];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nNLHHTalpha::ExternalNodeCondition(void) const
{
	return KBC_CardT::kAcc;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nNLHHTalpha::nComputeParameters(void)
{
	/* predictor */
	dpred_v = fdt;
	dpred_a	= (1.0 - 2.0*fbeta)*0.5*fdt*fdt;
	vpred_a	= (1.0 - fgamma)*fdt;
	
	/* corrector/consistent BC */
	dcorr_a = fbeta*fdt*fdt;
	vcorr_a = fgamma*fdt;
}
