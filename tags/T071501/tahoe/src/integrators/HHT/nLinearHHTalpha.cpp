/* $Id: nLinearHHTalpha.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "nLinearHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nLinearHHTalpha::nLinearHHTalpha(ifstreamT& in, ostream& out, int auto2ndorder):
	HHTalpha(in, out, auto2ndorder)
{

}

/* consistent BC's - updates predictors and acceleration only */
void nLinearHHTalpha::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nLinearHHTalpha::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU)(node, dof);
	double& v = (*fdU)(node, dof);
	double& a = (*fddU)(node, dof);

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
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nLinearHHTalpha::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU->Length() != fdU->Length() ||
	   fdU->Length() != fddU->Length()) throw eGeneralFail;
#endif
	
	/* save values from t_n (need by HHT-alpha) */
	dn = *fU;
	vn = *fdU;

	/* displacement predictor */
	fU->AddCombination(dpred_v, *fdU, dpred_a, *fddU);

	/* velocity predictor */
	fdU->AddScaled(vpred_a, *fddU);
}		

/* correctors - map ACTIVE */
void nLinearHHTalpha::Corrector(const iArray2DT& eqnos, const dArrayT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nLinearHHTalpha::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU->Length()   ||
	      fU->Length() != fdU->Length()  ||
	     fdU->Length() != fddU->Length() ||
		fddU->Length() != dn.Length()    ||
		   dn.Length() != vn.Length()) throw eGeneralFail;		
	/* note: no check on length of update.  could make numequations
	         a field in ControllerT. */
#endif

	/* displacement */
	(*fU) *= dcorr_dpred;
	fU->AddScaled(dcorr_d, dn);

	/* velocity */
	(*fdU) *= vcorr_vpred;
	fdU->AddScaled(vcorr_v, vn);

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = fU->Pointer();
	double *pv  = fdU->Pointer();
	double *pa  = fddU->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++;
		
		/* active dof */
		if (eq > 0)
		{
			double a = update[--eq]; //OFFSET
		
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
	const dArray2DT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nLinearHHTalpha::MappedCorrector: field arrays not set" << endl;
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
				double a = *pupdate;
			
				*pd += dcorr_a*a;
				*pv += vcorr_a*a;
				*pa = a;
			}
			
			/* next */
			pupdate++; pd++; pv++; pa++;
		}
	}
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
