/* $Id: nStaticIntegrator.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "nStaticIntegrator.h"
#include "ExceptionCodes.h"
#include "dArray2DT.h"

/* constructor */
nStaticIntegrator::nStaticIntegrator(void) { }

/* consistent BC's */
void nStaticIntegrator::ConsistentKBC(const KBC_CardT& KBC)
{
#if __option(extended_errorcheck)
	if (!fU)
	{
		cout << "\n nStaticIntegrator::ConsistentKBC: field arrays not set" << endl;
		throw eGeneralFail;
	}
#endif

	/* destination */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (*fU)(node, dof);
	
	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		{
			d = 0.0;
			break;
		}
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			d = KBC.Value();
			break;
		}
		default:
			cout << "\n nTrapezoid::ConsistentKBC:unknown BC code\n" << endl;
			throw eBadInputValue;
	}
}		

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nStaticIntegrator::ExternalNodeCondition(void) const
{
	return KBC_CardT::kDsp;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate time stepping constants */
void nStaticIntegrator::nComputeParameters(void) { }
