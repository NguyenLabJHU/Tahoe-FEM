/* $Id: PenaltySphereT.cpp,v 1.11 2004-07-15 08:31:15 paklein Exp $ */
/* created: paklein (04/30/1998) */
#include "PenaltySphereT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* constructor */
PenaltySphereT::PenaltySphereT(void):
	fLHS(ElementMatrixT::kSymmetric)
{
	SetName("sphere_penalty");
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
	
	/* equations */
	const iArray2DT& eqnos = Field().Equations();

	/* support class */
	const FieldSupportT& support = FieldSupport();

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
			eqnos.RowAlias(fContactNodes[i], fi_sh);
			support.AssembleLHS(fGroup, fLHS, fi_sh);
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

/* accept parameter list */
void PenaltySphereT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PenaltyRegionT::TakeParameterList(list);

	/* sphere dimensions */
	fRadius = list.GetParameter("radius");

	/* dimension work space */
	fv_OP.Dimension(Field().NumDOF());
	fLHS.Dimension(FieldSupport().NumSD());
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltySphereT::ComputeContactForce(double kforce)
{
	/* coordinates */
	const dArray2DT& coords = FieldSupport().CurrentCoordinates();

	/* loop over strikers */
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		coords.RowCopy(fContactNodes[i], fv_OP);
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
