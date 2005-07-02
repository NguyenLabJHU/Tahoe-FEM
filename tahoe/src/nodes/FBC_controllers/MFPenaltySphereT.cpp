/* $Id: MFPenaltySphereT.cpp,v 1.8.20.1 2005-07-02 22:50:38 paklein Exp $ */
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

	/* register with memory managers */
	fdContactNodesGroup2D.Register(fInitCoords);
	fdContactNodesGroup2D.Register(fCurrCoords);
}

/* (re-)set the configuration */
GlobalT::InitStatusT MFPenaltySphereT::UpdateConfiguration(void)
{
	/* inherited */
	GlobalT::InitStatusT status = PenaltySphereT::UpdateConfiguration();

	/* collect coordinates */
	fInitCoords.RowCollect(fContactNodes, FieldSupport().InitialCoordinates());

	return status;
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
	fCurrCoords += fInitCoords; //EFFECTIVE DVA

	/* loop over strikers */
	fContactForce2D = 0.0;	
	for (int i = 0; i < fContactNodes.Length(); i++)
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
