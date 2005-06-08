/* $Id: nNLHHTalpha.cpp,v 1.14.12.1 2005-06-08 17:22:50 paklein Exp $ */
/* created: paklein (10/17/1996) */
#include "nNLHHTalpha.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"
#include "BasicFieldT.h"

using namespace Tahoe;

/* constructor */
nNLHHTalpha::nNLHHTalpha(double alpha):
	HHTalpha(alpha)
{

}

/* consistent BC's - updates predictors and acceleration only */
void nNLHHTalpha::ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC, const iArrayT* nodes)
{
	const char caller[] = "nNLHHTalpha::ConsistentKBC";
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

		switch ( KBC.Code() )
		{
			case KBC_CardT::kFix: /* zero displacement */
			case KBC_CardT::kDsp: /* prescribed displacement */
			{
				double temp = KBC.Value();
				
				if (fabs(dcorr_a) > kSmall) /* for dt -> 0.0 */
					a = (temp - d)/dcorr_a;
				else
					a = 0.0;
				d  = temp;
				v += vcorr_a*a;

				break;
			}
			case KBC_CardT::kVel: /* prescribed velocity */
			{
				double temp = KBC.Value();
				
			if (fabs(vcorr_a) > kSmall) /* for dt -> 0.0 */
				a = (temp - v)/vcorr_a;
			else
				a = 0.0;
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
			case KBC_CardT::kNull: /* do nothing */
			{
				break;
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown BC code: %d", KBC.Code());
		}
	}
}
	
#pragma message ("roll up redundancy after it works")
// predictors - map ALL, unless limit arguments are specified
void nNLHHTalpha::Predictor(BasicFieldT& field, int fieldstart /*= 0*/, int fieldend /*= -1*/)
{
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
void nNLHHTalpha::Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart /*= 0*/, int fieldend /*= -1*/, int dummy /*= 0*/)
{
#pragma unused(dummy)

	if (fieldend == -1) // operate on full arrays
	{
		/* displacement corrector */
		field[0].AddScaled(dcorr_a, update);
		
		/* velocity corrector */
		field[1].AddScaled(vcorr_a, update);
		
		/* acceleration corrector */
		field[2] += update;
	}
	else // operate on restricted contiguous block of the arrays
	{
		/* displacement corrector */
		field[0].AddScaled(dcorr_a, update, fieldstart, fieldend);
		
		/* velocity corrector */
		field[1].AddScaled(vcorr_a, update, fieldstart, fieldend);
		
		/* acceleration corrector */
		field[2].AddScaled(1.0, update, fieldstart, fieldend);
	}
}

/* correctors - map ACTIVE */
void nNLHHTalpha::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	const iArray2DT& eqnos = field.Equations();

	/* add update - assumes that fEqnos maps directly into dva */
	const int *peq = eqnos.Pointer();
	
	double *pd  = field[0].Pointer();
	double *pv  = field[1].Pointer();
	double *pa  = field[2].Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++ - eq_start;
		
		if (eq > -1 && eq < num_eq) /* active dof */
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

void nNLHHTalpha::MappedCorrector(BasicFieldT& field, const iArrayT& map, 
		const iArray2DT& flags, const dArray2DT& update)
{
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
				double da = *pupdate - *pa; 
/* NOTE: update is the total a_n+1, so we need to recover da, the acceleration
 *       increment. Due to truncation errors, this will not match the update
 *       applied in nNLHHTalpha::Corrector exactly. Therefore, ghosted nodes
 *       will not have exactly the same trajectory as their images. The solution
 *       would be to expand da to the full field, and send it. Otherwise, we
 *       could develop maps from the update vector to each communicated outgoing
 *       packet. Or {d,v,a} could be recalculated from t_n here and in Corrector,
 *       though this would require the data from t_n. */
			
				*pd += dcorr_a*da;
				*pv += vcorr_a*da;
				*pa  = *pupdate; /* a_n+1 matches exactly */
			}
			
			/* next */
			pflag++; pupdate++; pd++; pv++; pa++;
		}
	}
}

/* return the field array needed by nIntegratorT::MappedCorrector. */
const dArray2DT& nNLHHTalpha::MappedCorrectorField(BasicFieldT& field) const
{
	return field[2];
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
