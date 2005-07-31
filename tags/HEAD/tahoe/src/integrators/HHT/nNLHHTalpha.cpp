/* $Id: nNLHHTalpha.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/17/1996)                                          */

#include "nNLHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nNLHHTalpha::nNLHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder):
	HHTalpha(in, out, auto2ndorder)
{

}

/* consistent BC's - updates predictors and acceleration only */
void nNLHHTalpha::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nNLHHTalpha::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU)(node, dof);
	double& v = (*fdU)(node, dof);
	double& a = (*fddU)(node, dof);

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
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nNLHHTalpha::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU->Length() != fdU->Length() ||
	   fdU->Length() != fddU->Length()) throw eGeneralFail;
#endif

	/* displacement predictor */
	fU->AddCombination(dpred_v, *fdU, dpred_a, *fddU);

	/* velocity predictor */
	fdU->AddScaled(vpred_a, *fddU);
	
	/* acceleration predictor */
	(*fddU) = 0.0;
}		

/* correctors - map ACTIVE */
void nNLHHTalpha::Corrector(const iArray2DT& eqnos, const dArrayT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nNLHHTalpha::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU->Length()  ||
	      fU->Length() != fdU->Length() ||
	     fdU->Length() != fddU->Length()) throw eGeneralFail;		
	/* note: No check on length of update.  Could make numequations
	         a field in ControllerT. */
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	
	double *pd  = fU->Pointer();
	double *pv  = fdU->Pointer();
	double *pa  = fddU->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++;
		
		if (eq > 0)	/* active dof */
		{
			double da = update[--eq];
		
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
	const dArray2DT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nNLHHTalpha::MappedCorrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* checks */
	if (flags.MajorDim() != map.Length() ||
	    flags.MajorDim() != update.MajorDim() ||
	    flags.MinorDim() != update.MinorDim() ||
	    flags.MinorDim() != fU->MinorDim()) throw eSizeMismatch;

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		int* pflags = flags(i);

		double* pd = (*fU)(row);
		double* pv = (*fdU)(row);
		double* pa = (*fddU)(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag++ > 0)
			{
				double da = *pupdate;
			
				*pd += dcorr_a*da;
				*pv += vcorr_a*da;
				*pa += da;
			}
			
			/* next */
			pupdate++; pd++; pv++; pa++;
		}
	}
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
