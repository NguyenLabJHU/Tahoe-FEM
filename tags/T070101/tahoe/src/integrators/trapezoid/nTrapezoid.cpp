/* $Id: nTrapezoid.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/03/1999)                                          */

#include "nTrapezoid.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"

/* constructor */
nTrapezoid::nTrapezoid(void) { }

/* consistent BC's */
void nTrapezoid::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU || !fdU)
	{
		cout << "\n nTrapezoid::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU)(node, dof);
	double& v = (*fdU)(node, dof);
	
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
	if (!fU || !fdU)
	{
		cout << "\n nTrapezoid::Predictor: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (fU->Length() != fdU->Length()) throw eSizeMismatch;
#endif
	
	/* displacement predictor */
	fU->AddScaled(dpred_v, *fdU);
}		

/* correctors - map ACTIVE */
void nTrapezoid::Corrector(const iArray2DT& eqnos, const dArrayT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU)
	{
		cout << "\n nTrapezoid::Corrector: field arrays not set" << endl;
		throw eGeneralFail;
	}
	if (eqnos.Length() != fU->Length()   ||
	      fU->Length() != fdU->Length()) throw eSizeMismatch;		
	//NOTE: no check on length of update.  could make numequations
	//      a field in ControllerT.
#endif

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	double *pd  = fU->Pointer();
	double *pv  = fdU->Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = *peq++;
		/* active dof */
		if (eq > 0)
		{
			double v = update[--eq]; //OFFSET
			*pd += dcorr_v*v;
			*pv += v;
		}
		pd++;
		pv++;
	}
}

void nTrapezoid::MappedCorrector(const iArrayT& map, const iArray2DT& flags,
	const dArray2DT& update)
{
#if __option (extended_errorcheck)
	if (!fU || !fdU)
	{
		cout << "\n nTrapezoid::MappedCorrector: field arrays not set" << endl;
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
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag++ > 0)
			{
				double v = *pupdate;
				*pd += dcorr_v*v;
				*pv += v;
			}
			
			/* next */
			pupdate++; pd++; pv++;
		}
	}
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
