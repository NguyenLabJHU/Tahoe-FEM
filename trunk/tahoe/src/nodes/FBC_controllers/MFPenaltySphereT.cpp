/* $Id: MFPenaltySphereT.cpp,v 1.8 2004-07-15 08:31:15 paklein Exp $ */
/* created: paklein (04/17/2000) */
#include "MFPenaltySphereT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "ElementBaseT.h"

using namespace Tahoe;

/* constructor */
MFPenaltySphereT::MFPenaltySphereT(void):
	fElementGroup(NULL)
{
	SetName("sphere_penalty_meshfree");
}

/* describe the parameters needed by the interface */
void MFPenaltySphereT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PenaltySphereT::DefineParameters(list);

	/* meshless group number */
	list.AddParameter(ParameterT::Integer, "meshless_element_group");
}

/* accept parameter list */
void MFPenaltySphereT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFPenaltySphereT::TakeParameterList";

	/* inherited */
	PenaltySphereT::TakeParameterList(list);

	/* meshless group number */
	int group = list.GetParameter("meshless_element_group");
	group--;
	/* get pointer */
	fElementGroup = &(FieldSupport().ElementGroup(group));
	if (!fElementGroup)
		ExceptionT::GeneralFail(caller, "error retrieving pointer for group %d", group+1);
	
	/* check */
	if (fElementGroup->InterpolantDOFs())
		ExceptionT::GeneralFail(caller, "element group %d has interpolant DOF's. Use PenaltySphereT", group+1);

	/* allocate workspace */
	const dArray2DT& coords = FieldSupport().CurrentCoordinates();
	fCoords.Dimension(fNumContactNodes, coords.MinorDim());
	fCurrCoords.Dimension(fNumContactNodes, Field().NumDOF());
	
	/* collect coordinates */
	fCoords.RowCollect(fContactNodes, coords);
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* compute the nodal contribution to the residual force vector */
void MFPenaltySphereT::ComputeContactForce(double kforce)
{
	/* DOF source */
	if (!fElementGroup) SetElementGroup();

	/* get nodal displacements */
	fElementGroup->NodalDOFs(fContactNodes, fCurrCoords);
	fCurrCoords += fCoords; //EFFECTIVE DVA

	/* loop over strikers */
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		fCurrCoords.RowCopy(i, fv_OP);
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

/**********************************************************************
* Private
**********************************************************************/

/* get element group pointer */
void MFPenaltySphereT::SetElementGroup(void)
{
}
