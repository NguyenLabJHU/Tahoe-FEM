/* $Id: nExplicitCD.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Node controller for an explicit 2nd order                              */
/* accurate, central difference time-stepping algorithm.                  */

#include "nExplicitCD.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nExplicitCD::nExplicitCD(void) { }

/* consistent BC's - updates predictors and acceleration only */
void nExplicitCD::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nExplicitCD::ConsistentKBC: field arrays not set" << endl;
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
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nExplicitCD::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU->Length() != fdU->Length() ||
	   fdU->Length() != fddU->Length()) throw eGeneralFail;
#endif

	/* displacement predictor */
	fU->AddCombination(dpred_v, *fdU, dpred_a, *fddU);	

	/* velocity predictor */
	fdU->AddScaled(vpred_a, *fddU);
}		

/* correctors - map ACTIVE */
void nExplicitCD::Corrector(const iArray2DT& eqnos, const dArrayT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nExplicitCD::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU->Length() ||
	      fU->Length() != fdU->Length()  ||
	     fdU->Length() != fddU->Length()) throw eGeneralFail;		
	//NOTE: no check on length of update.  could make numequations
	//      a field in ControllerT.
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pv  = fdU->Pointer();
	double *pa  = fddU->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++;
		
		/* active dof */
		if (eq > 0)
		{
			double a = update[--eq]; //OFFSET
			*pv += vcorr_a*a;
			*pa = a;
		}
		pv++;
		pa++;
	}
}

void nExplicitCD::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU || !fddU)
	{
		cout << "\n nExplicitCD::MappedCorrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* checks */
	if (flags.MajorDim() != map.Length() ||
	    flags.MajorDim() != update.MajorDim() ||
	    flags.MinorDim() != update.MinorDim() ||
	    flags.MinorDim() != fdU->MinorDim()) throw eSizeMismatch;

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		int* pflags = flags(i);

		double* pv = (*fdU)(row);
		double* pa = (*fddU)(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag++ > 0)
			{
				double a = *pupdate;
				*pv += vcorr_a*a;
				*pa = a;
			}
			
			/* next */
			pupdate++; pv++; pa++;
		}
	}
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
