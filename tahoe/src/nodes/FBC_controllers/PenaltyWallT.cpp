/* $Id: PenaltyWallT.cpp,v 1.15 2005-11-07 21:00:23 regueiro Exp $ */
/* created: paklein (02/25/1997) */
#include "PenaltyWallT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "Vector3T.h"

using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
PenaltyWallT::PenaltyWallT(void):
	fLHS(ElementMatrixT::kSymmetric)
{
	SetName("wall_penalty");
}

/* tangent */
void PenaltyWallT::ApplyLHS(GlobalT::SystemTypeT sys_type)
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
		/* retrieve normal component */
		double normal_comp = fp_i.DotRow(i, fnormal);

		/* active */
		if (normal_comp < 0.0)
		{
			/* get force vector (has normal direction) */
			fContactForce2D.RowAlias(i, fd_sh);
			fLHS.Outer(fd_sh, fd_sh);
			fLHS *= constK/(normal_comp*normal_comp*fk);		
					
			/* assemble */
			eqnos.RowAlias(fContactNodes[i], fi_sh);
			support.AssembleLHS(fGroup, fLHS, fi_sh);
		}
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltyWallT::ComputeContactForce(double kforce)
{
	const char caller[] = "PenaltyWallT::ComputeContactForce";

	/* with "friction */
	if (fmu > kSmall)
	{
		//TEMP
		//ExceptionT::GeneralFail(caller, "general (2D/3D) friction implementation is not available");
//#if 0
		/* compute relative positions and velocities */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		fp_i.RowCollect(fContactNodes, coords);
		//where are velocities calculated?
		const dArray2DT& vels = FieldSupport().CurrentCoordinates();
		fv_i.RowCollect(fContactNodes, vels);
		for (int j = 0; j < fNumContactNodes; j++)
		{
			fp_i.AddToRowScaled(j, -1.0, fx);
			fv_i.AddToRowScaled(j, -1.0, fv);
		}
	
		/* compute contact forces */
		fntforce = 0.0;
		fContactForce2D = 0.0;	
		for (int i = 0; i < fNumContactNodes; i++)
		{
			double normal_comp = fp_i.DotRow(i, fnormal);
		
			/* penetration */
			if (normal_comp < 0.0)
			{			
				/* calculate tangential force due to friction */
				double tangent_comp = fv_i.DotRow(i, ftangent);
				fntforce[0] =-fk*normal_comp*kforce;		
				fntforce[1] = ((tangent_comp > 0.0) ? -1.0 : 1.0)*fmu*fntforce[0]*kforce;

				/* transform to x-y coordinates */
				fQ.Multx(fntforce, fxyforce);
				
				fContactForce2D.SetRow(i,fxyforce);
			}
			
			/* store gap */
			fGap[i] = normal_comp;
		}
//#endif		
	}
	else
	{
		/* compute relative positions */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		fp_i.RowCollect(fContactNodes, coords);
		for (int j = 0; j < fNumContactNodes; j++)
			fp_i.AddToRowScaled(j, -1.0, fx);
	
		/* compute contact forces */
		fntforce = 0.0;
		fContactForce2D = 0.0;	
		for (int i = 0; i < fNumContactNodes; i++)
		{
			double normal_comp = fp_i.DotRow(i, fnormal);
		
			/* penetration */
			if (normal_comp < 0.0)
			{
				/* normal force */
				fntforce[0] =-fk*normal_comp*kforce;		
		
				/* transform to x-y coordinates */
				fQ.Multx(fntforce, fxyforce);
				
				fContactForce2D.SetRow(i, fxyforce);
			}

			/* store gap */
			fGap[i] = normal_comp;
		}
	}
}

/* information about subordinate parameter lists */
void PenaltyWallT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	PenaltyRegionT::DefineSubs(sub_list);
	
	/* normal to the wall */
	sub_list.AddSub("wall_normal");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PenaltyWallT::NewSub(const StringT& name) const
{
	if (name == "wall_normal") {

		ParameterContainerT* n_choice = new ParameterContainerT(name);
		
		/* by dimension */
		n_choice->SetListOrder(ParameterListT::Choice);
		n_choice->AddSub("Vector_2");
		n_choice->AddSub("Vector_3");

		return n_choice;	
	}
	else /* inherited */
		return PenaltyRegionT::NewSub(name);
}

/* accept parameter list */
void PenaltyWallT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PenaltyWallT::TakeParameterList";

	/* inherited */
	PenaltyRegionT::TakeParameterList(list);

	/* get normal */
	int nsd = FieldSupport().NumSD();
	const ParameterListT& normal = list.GetListChoice(*this, "wall_normal");
	VectorParameterT::Extract(normal, fnormal);
	fnormal.UnitVector();
	if (fnormal.Length() != nsd) 
		ExceptionT::GeneralFail(caller, "\"wall_normal\" should be length %d not %d", nsd, fnormal.Length());
	
	/* transformation tensor */
	fQ.Dimension(nsd);
	if (nsd == 2)
	{
		/* set column vectors */
		fQ.SetCol(0, fnormal);
		fQ(0,1) =-fnormal[1];
		fQ(1,1) = fnormal[0];
	}
	else if (nsd == 3)
	{
		Vector3T<double> v0(fnormal.Pointer()), v1, v2; 
		
		/* find non-colinear directions */
		int i = 0;
		v2.Random(++i);
		while (v2.Norm() < 1.0e-06 || Vector3T<double>::Dot(v0,v2) > 0.99)
			v2.Random(++i);	
		v1.Cross(v0,v2);
		v2.Cross(v0,v1);

		/* write column vectors */
		fQ.SetCol(0,v0);
		fQ.SetCol(1,v1);
		fQ.SetCol(2,v2);
	}
	else 
		ExceptionT::GeneralFail(caller);
	
	/* dimension work space */	
	fntforce.Dimension(nsd);
	fxyforce.Dimension(nsd);
	fp_i.Dimension(fNumContactNodes, nsd);
	fv_i.Dimension(fNumContactNodes, nsd);
	fLHS.Dimension(nsd);
}
