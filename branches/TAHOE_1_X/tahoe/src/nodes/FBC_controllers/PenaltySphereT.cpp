/* $Id: PenaltySphereT.cpp,v 1.10 2004-06-17 07:41:53 paklein Exp $ */
/* created: paklein (04/30/1998) */
#include "PenaltySphereT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "ifstreamT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* constructor */
PenaltySphereT::PenaltySphereT(FEManagerT& fe_manager,
	int group,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT& disp,
	const dArray2DT* vels):
	PenaltyRegionT(fe_manager, group, eqnos, coords, disp, vels),
	fv_OP(rCoords.MinorDim()),
	fLHS(eqnos.MinorDim(),ElementMatrixT::kSymmetric)
{
	SetName("sphere_penalty");
}

/* input processing */
void PenaltySphereT::EchoData(ifstreamT& in, ostream& out)
{
	/* inherited */
	PenaltyRegionT::EchoData(in, out);

	/* echo parameters */
	in >> fRadius; if (fRadius < 0.0) throw ExceptionT::kBadInputValue;

	out << " Sphere radius . . . . . . . . . . . . . . . . . = " << fRadius << '\n';
}

/* form of tangent matrix */
GlobalT::SystemTypeT PenaltySphereT::TangentType(void) const
{
	/* symmetric tangent (frictionless) */
	return GlobalT::kSymmetric;
}

/* tangent term */
void PenaltySphereT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* node by node */
	for (int i = 0; i < fNumContactNodes; i++)
	{
		double gap  = fGap[i];
		double dist = gap + fRadius;

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
			fFEManager.AssembleLHS(fGroup, fLHS, fi_sh);
		}
	}
}

/* describe the parameters needed by the interface */
void PenaltySphereT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PenaltyRegionT::DefineParameters(list);
	
	list.AddParameter(fRadius, "radius");
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltySphereT::ComputeContactForce(double kforce)
{
	/* loop over strikers */
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		rCoords.RowCopy(fContactNodes[i], fv_OP);
		fv_OP -= fx;
		
		/* penetration */
		double dist = fv_OP.Magnitude();
		double pen  = dist - fRadius;
		if (pen < 0.0)
		{
			/* convert to force*outward normal */
			fv_OP *= (-pen*fk*kforce/dist);
		
			/* accumulate */
			fContactForce2D.SetRow(i, fv_OP);
		}

		/* store */
		fGap[i] = pen;
	}
}
