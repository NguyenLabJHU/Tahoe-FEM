/* $Id: PenaltyCylinderT.cpp,v 1.4 2004-07-15 08:31:15 paklein Exp $ */
#include "PenaltyCylinderT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* constructor */
PenaltyCylinderT::PenaltyCylinderT(void):
	fLHS(ElementMatrixT::kSymmetric)
{
	SetName("cylinder_penalty");
}

/* form of tangent matrix */
GlobalT::SystemTypeT PenaltyCylinderT::TangentType(void) const
{
	/* symmetric tangent (frictionless) */
	return GlobalT::kSymmetric;
}

/* tangent term */
void PenaltyCylinderT::ApplyLHS(GlobalT::SystemTypeT sys_type)
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
			fLHS.Outer(fd_sh, fd_sh, constK*((fk/dPhi) - (1.0/dist))/dPhi);
			fLHS.PlusIdentity(constK*dPhi/dist);
			fLHS.Outer(fDirection, fDirection, -constK*dPhi/dist, dMatrixT::kAccumulate);
		
			/* assemble */
			eqnos.RowAlias(fContactNodes[i], fi_sh);
			support.AssembleLHS(fGroup, fLHS, fi_sh);
		}
	}
}

/* describe the parameters needed by the interface */
void PenaltyCylinderT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PenaltyRegionT::DefineParameters(list);
	
	list.AddParameter(fRadius, "radius");
}

/* information about subordinate parameter lists */
void PenaltyCylinderT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	PenaltyRegionT::DefineSubs(sub_list);

	/* direction */
	sub_list.AddSub("cylinder_penalty_axis");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PenaltyCylinderT::NewSub(const StringT& name) const
{
	if (name == "cylinder_penalty_axis")
	{
		ParameterContainerT* dir = new ParameterContainerT(name);
		
		/* by dimension */
		dir->SetListOrder(ParameterListT::Choice);
		dir->AddSub("Vector_2");
		dir->AddSub("Vector_3");
	
		return dir;
	}
	else /* inherited */
		return PenaltyRegionT::NewSub(name);
}

/* accept parameter list */
void PenaltyCylinderT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PenaltyCylinderT::TakeParameterList";

	/* inherited */
	PenaltyRegionT::TakeParameterList(list);

	/* get parameters */
	fRadius = list.GetParameter("radius");

	/* initial position */
	int nsd = FieldSupport().NumSD();
	const ParameterListT& dir = list.GetListChoice(*this, "cylinder_penalty_axis");
	VectorParameterT::Extract(dir, fDirection);
	fDirection.UnitVector();
	if (fDirection.Length() != nsd) 
		ExceptionT::GeneralFail(caller, "\"cylinder_penalty_axis\" should be length %d not %d", nsd, fDirection.Length());

	/* dimension work space */
	fDirection.Dimension(nsd);
	fR.Dimension(nsd);
	fv_OP.Dimension(nsd);
	fLHS.Dimension(nsd);
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltyCylinderT::ComputeContactForce(double kforce)
{
	/* loop over strikers */
	const dArray2DT& coords = FieldSupport().CurrentCoordinates();
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		coords.RowCopy(fContactNodes[i], fv_OP);
		fv_OP -= fx;

		/* vector in radial direction */
		fR.SetToCombination(1.0, fv_OP, -dArrayT::Dot(fv_OP, fDirection), fDirection);
		
		/* penetration */
		double dist = fR.Magnitude();
		double pen  = dist - fRadius;
		if (pen < 0.0)
		{
			/* convert to force*outward normal */
			fR *= (-pen*fk*kforce/dist);
		
			/* accumulate */
			fContactForce2D.SetRow(i, fR);
		}

		/* store */
		fGap[i] = pen;
	}
}
