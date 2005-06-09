/* $Id: nLinearHHTalpha.cpp,v 1.14.12.2 2005-06-09 02:43:53 d-farrell2 Exp $ */
/* created: paklein (10/14/1996) */
#include "nLinearHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"
#include "BasicFieldT.h"

using namespace Tahoe;

/* constructor */
nLinearHHTalpha::nLinearHHTalpha(double alpha):
	HHTalpha(alpha),
	fField(NULL)
{

}

/* register field with the integrator */
void nLinearHHTalpha::Dimension(const BasicFieldT& field)
{
	/* inherited */
	nIntegratorT::Dimension(field);

	//TEMP - can only handle single field
	if (!fField) 
		fField = &field;
	else if (fField != &field)
		ExceptionT::GeneralFail("nLinearHHTalpha::Dimension", "only single field allowed");

	/* dimension saved values from t_n (need by HHT-alpha) */
	dn.Dimension(field[0]);
	dn = 0.0;
	vn.Dimension(field[1]);
	vn = 0.0;
}

/* prescribe the field and derivatives consistent BC's */
void nLinearHHTalpha::ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC, const iArrayT* nodes)
{
	const char caller[] = "nLinearHHTalpha::ConsistentKBC";
	if (KBC.Mode() == KBC_CardT::kSet && !nodes)
		ExceptionT::GeneralFail(caller, "expecting non-NULL nodes");

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	int nnd = (KBC.Mode() == KBC_CardT::kSet) ? nodes->Length() : 1;
	const int* pnd = (KBC.Mode() == KBC_CardT::kSet) ? nodes->Pointer() : &node;

	/* apply to nodes */
	for (int i = 0; i < nnd; i++)
	{
		double& d = (field[0])(*pnd, dof);
		double& v = (field[1])(*pnd, dof);
		double& a = (field[2])(*pnd, dof);
		pnd++; /* next */

		switch (KBC.Code())
		{
			case KBC_CardT::kFix: /* zero displacement */
			case KBC_CardT::kDsp: /* prescribed displacement */
			{
				double d_alpha = (1.0 + falpha)*KBC.Value() - falpha*dn(node,dof);

				if (fabs(dalpha_a) > kSmall) /* for dt -> 0.0 */
					a = (d_alpha - d)/dalpha_a;
				else
					a = 0.0;
				d  = d_alpha;
				v += valpha_a*a;
		
				break;
			}
			case KBC_CardT::kVel: /* prescribed velocity */
			{
				double v_alpha = (1.0 + falpha)*KBC.Value() - falpha*vn(node,dof);
	
				if (fabs(valpha_a) > kSmall) /* for dt -> 0.0 */
					a = (v_alpha - v)/valpha_a;
				else
					a = 0.0;
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
			case KBC_CardT::kNull: /* do nothing */
				break;
			default:
				ExceptionT::BadInputValue(caller, "unknown BC code: %d", KBC.Code());
		}		
	}
}		
#pragma message ("roll up redundancy after it works")
// predictors - map ALL, unless limit arguments are specified
void nLinearHHTalpha::Predictor(BasicFieldT& field, int fieldstart /*= 0*/, int fieldend /*= -1*/)
{
	/* save values from t_n (need by HHT-alpha) */
	dn = field[0];
	vn = field[1];

	if (fieldend == -1) // operate on full arrays
	{
		/* displacement predictor */
		field[0].AddCombination(dpred_v, field[1], dpred_a, field[2]);
		
		/* velocity predictor */
		field[1].AddScaled(vpred_a, field[2]);
		
		/* acceleratior predictor */
		field[2] = 0.0;	
	}
	else // operate on restricted contiguous block of the arrays
	{
		/* displacement predictor */
		field[0].AddCombination(dpred_v, field[1], dpred_a, field[2], fieldstart, fieldend);
		
		/* velocity predictor */
		field[1].AddScaled(vpred_a, field[2], fieldstart, fieldend);
		
		/* acceleratior predictor */
		field[2].SetToScaled(0.0, field[1], fieldstart, fieldend);	
	}
}		

/* corrector. Maps ALL degrees of freedom forward. */
void nLinearHHTalpha::Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart /*= 0*/, int fieldend /*= -1*/)
{

	if (fieldend == -1) // operate on full arrays
	{
		/* displacement corrector */
		field[0] *= dcorr_dpred;
		field[0].AddCombination(dcorr_d, dn, dcorr_a, update);
		
		/* velocity corrector */
		field[1] *= vcorr_vpred;
		field[1].AddCombination(vcorr_v, vn, vcorr_a, update);
		
		/* acceleration corrector */
		field[2] = update;
	}
	else // operate on restricted contiguous block of the arrays
	{
		/* displacement corrector */
		field[0].SetToScaled(dcorr_dpred, field[0], fieldstart, fieldend);
		field[0].AddCombination(dcorr_d, dn, dcorr_a, update, fieldstart, fieldend);
		
		/* velocity corrector */
		field[1].SetToScaled(vcorr_vpred, field[1], fieldstart, fieldend);
		field[1].AddCombination(vcorr_v, vn, vcorr_a, update, fieldstart, fieldend);
		
		/* acceleration corrector */
		field[2].SetToScaled(1.0, update, fieldstart, fieldend);
	}
#pragma message("Not exctly sure about this one")
}

/* correctors - map ACTIVE */
void nLinearHHTalpha::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	/* displacement */
	field[0] *= dcorr_dpred;
	field[0].AddScaled(dcorr_d, dn);

	/* velocity */
	field[1] *= vcorr_vpred;
	field[1].AddScaled(vcorr_v, vn);

	/* add update - assumes that fEqnos maps directly into dva */
	const iArray2DT& eqnos = field.Equations();
	
	const int *peq = eqnos.Pointer();
	double *pd = field[0].Pointer();
	double *pv = field[1].Pointer();
	double *pa = field[2].Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		
		/* active dof */
		if (eq > -1 && eq < num_eq)
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

void nLinearHHTalpha::MappedCorrector(BasicFieldT& field, const iArrayT& map,
		const iArray2DT& flags, const dArray2DT& update)
{
	/* check dimensions */
	if (flags.MajorDim() != update.MajorDim() ||
	    flags.MinorDim() != update.MinorDim()) throw ExceptionT::kSizeMismatch;

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		const int* pflags = flags(i);

		double* pd = (field[0])(row);
		double* pv = (field[1])(row);
		double* pa = (field[2])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag > 0)
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

/* return the field array needed by nIntegratorT::MappedCorrector. */
const dArray2DT& nLinearHHTalpha::MappedCorrectorField(BasicFieldT& field) const
{
	return field[2];
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
