/* $Id: PenaltySphereT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/30/1998)                                          */

#include "PenaltySphereT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "FEManagerT.h"
#include "fstreamT.h"
#include "eControllerT.h"

/* constructor */
PenaltySphereT::PenaltySphereT(FEManagerT& fe_manager,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	PenaltyRegionT(fe_manager, eqnos, coords, vels),
	fv_OP(rCoords.MinorDim()),
	fLHS(eqnos.MinorDim(),ElementMatrixT::kSymmetric)
{

}

/* input processing */
void PenaltySphereT::EchoData(ifstreamT& in, ostream& out)
{
	/* inherited */
	PenaltyRegionT::EchoData(in, out);

	/* echo parameters */
	in >> fRadius; if (fRadius < 0.0) throw eBadInputValue;

	out << " Sphere radius . . . . . . . . . . . . . . . . . = " << fRadius << '\n';
}

/* initialize data */
void PenaltySphereT::Initialize(void)
{
	/* inherited */
	PenaltyRegionT::Initialize();
	
	/* memory for distances */
	fDistances.Allocate(fNumContactNodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT PenaltySphereT::TangentType(void) const
{
	/* symmetric tangent (frictionless) */
	return GlobalT::kSymmetric;
}

/* tangent term */
void PenaltySphereT::ApplyLHS(void)
{
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

	/* node by node */
	for (int i = 0; i < fNumContactNodes; i++)
	{
		double dist = fDistances[i];
		double gap  = dist - fRadius;

		/* active */
		if (gap < 0.0)
		{
			/* get force vector (has normal direction) */
			fContactForce2D.RowAlias(i,fd_sh);
			
			double dPhi = gap*fk;
			fLHS.Outer(fd_sh,fd_sh);
			fLHS *= constK*((fk/dPhi) - (1.0/dist))/dPhi;
		
			fLHS.PlusIdentity(constK*dPhi/dist);
		
			/* assemble */
			rEqnos.RowAlias(fContactNodes[i], fi_sh);
			fFEManager.AssembleLHS(fLHS, fi_sh);
		}
	}
}

/**********************************************************************
* Protected
**********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltySphereT::ComputeContactForce(double kforce)
{
	/* loop over strikers */
	fh_max = 0.0;
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		rCoords.RowCopy(fContactNodes[i], fv_OP);
		fv_OP -= fx;
		
		/* penetration */
		double dist = fv_OP.Magnitude();
		double pen  = fRadius - dist;
		if (pen > 0.0)
		{
			/* store max penetration */
			fh_max = (pen > fh_max) ? pen : fh_max;
		
			/* convert to force*outward normal */
			fv_OP *= (pen*fk*kforce/dist);
		
			/* accumulate */
			fContactForce2D.SetRow(i, fv_OP);
		}

		/* store */
		fDistances[i] = dist;
	}
}
