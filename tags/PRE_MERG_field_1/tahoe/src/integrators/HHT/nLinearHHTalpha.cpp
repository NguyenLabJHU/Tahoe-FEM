/* $Id: nLinearHHTalpha.cpp,v 1.3 2002-04-02 23:19:20 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "nLinearHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nLinearHHTalpha::nLinearHHTalpha(ifstreamT& in, ostream& out, bool auto2ndorder):
	HHTalpha(in, out, auto2ndorder),
	nControllerT(2)
{

}

/* consistent BC's - updates predictors and acceleration only */
void nLinearHHTalpha::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nLinearHHTalpha::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU[0])(node, dof);
	double& v = (*fU[1])(node, dof);
	double& a = (*fU[2])(node, dof);

	switch (KBC.Code())
	{
		case KBC_CardT::kFix: /* zero displacement */
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			double d_alpha = (1.0 + falpha)*KBC.Value() - falpha*dn(node,dof);

			a  = (d_alpha - d)/dalpha_a;
			d  = d_alpha;
			v += valpha_a*a;
		
			break;
		}
		
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			double v_alpha = (1.0 + falpha)*KBC.Value() - falpha*vn(node,dof);
	
			a  = (v_alpha - v)/valpha_a;
			v  = v_alpha;
			d += dalpha_a*a;
	
			break;
		}
		
		case KBC_CardT::kAcc: /* prescribed acceleration */
		{
			a  = KBC.Value();
			d += dalpha_a*a;
			v += valpha_a*a;

			break;
		}

		default:
		
			cout << "\nnLinearHHTalpha::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* predictors - map ALL */
void nLinearHHTalpha::Predictor(void)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nLinearHHTalpha::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU[0]->Length() != fU[1]->Length() ||
	    fU[1]->Length() != fU[2]->Length()) throw eGeneralFail;
#endif
	
	/* save values from t_n (need by HHT-alpha) */
	dn = *fU[0];
	vn = *fU[1];

	/* displacement predictor */
	fU[0]->AddCombination(dpred_v, *fU[1], dpred_a, *fU[2]);

	/* velocity predictor */
	fU[1]->AddScaled(vpred_a, *fU[2]);
}		

/* correctors - map ACTIVE */
void nLinearHHTalpha::Corrector(const iArray2DT& eqnos, const dArrayT& update,
	int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nLinearHHTalpha::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU[0]->Length()   ||
	   fU[0]->Length() != fU[1]->Length()  ||
	   fU[1]->Length() != fU[2]->Length() ||
       fU[2]->Length() != dn.Length()    ||
		   dn.Length() != vn.Length()) throw eGeneralFail;		
	/* note: no check on length of update. */
#endif

	/* displacement */
	(*fU[0]) *= dcorr_dpred;
	fU[0]->AddScaled(dcorr_d, dn);

	/* velocity */
	(*fU[1]) *= vcorr_vpred;
	fU[1]->AddScaled(vcorr_v, vn);

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = fU[0]->Pointer();
	double *pv  = fU[1]->Pointer();
	double *pa  = fU[2]->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		
		/* active dof */
		if (eq > -1 && eq < eq_stop)
		{
			double a = update[eq];
		
			*pd += dcorr_a*a;
			*pv += vcorr_a*a;
			*pa = a;
		}
			
		pd++;
		pv++;
		pa++;
	}

}

void nLinearHHTalpha::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update, int eq_start, int eq_stop)
{
#if __option (extended_errorcheck)
	if (fU[0] == NULL ||
	    fU[1] == NULL ||
	    fU[2] == NULL)
	{
		cout << "\n nLinearHHTalpha::MappedCorrector: field arrays not set" << endl;
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
				double a = *pupdate;
			
				*pd += dcorr_a*a;
				*pv += vcorr_a*a;
				*pa = a;
			}
			
			/* next */
			pflag++; pupdate++; pd++; pv++; pa++;
		}
	}
}

/* return the field array needed by nControllerT::MappedCorrector. */
const dArray2DT& nLinearHHTalpha::MappedCorrectorField(void) const
{
#if __option (extended_errorcheck)
	if (fU[2] == NULL)
	{
		cout << "\n nLinearHHTalpha::MappedCorrectorField: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	return *fU[2];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nLinearHHTalpha::ExternalNodeCondition(void) const
{
	return KBC_CardT::kAcc;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nLinearHHTalpha::nComputeParameters(void)
{
	/* predictor */
	dpred_v		= (1.0 + falpha)*fdt;
	dpred_a		= (1.0 + falpha)*(1.0 - 2.0*fbeta)*0.5*fdt*fdt;
	vpred_a		= (1.0 + falpha)*(1.0 - fgamma)*fdt;
	
	/* corrector */
	dcorr_d		= falpha/(1.0 + falpha);
	dcorr_dpred = 1.0/(1.0 + falpha);
	dcorr_a		= fbeta*fdt*fdt;
	
	vcorr_v		= dcorr_d;
	vcorr_vpred = dcorr_dpred;
	vcorr_a		= fgamma*fdt;

	/* consistent BC */
	dalpha_a    = (1.0 + falpha)*fbeta*fdt*fdt;
	valpha_a    = (1.0 + falpha)*fgamma*fdt;
}
