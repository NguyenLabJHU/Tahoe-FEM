/* $Id: MFPenaltySphereT.cpp,v 1.7 2004-06-17 07:41:53 paklein Exp $ */
/* created: paklein (04/17/2000) */
#include "MFPenaltySphereT.h"

#include <iostream.h>

#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "ElementBaseT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
MFPenaltySphereT::MFPenaltySphereT(FEManagerT& fe_manager, int group,
	const iArray2DT& eqnos, const dArray2DT& coords, const dArray2DT& disp, const dArray2DT* vels):
	PenaltySphereT(fe_manager, group, eqnos, coords, disp, vels),
	fElementGroup(NULL)
{

}

/* input processing */
void MFPenaltySphereT::EchoData(ifstreamT& in, ostream& out)
{
	/* inherited */
	PenaltySphereT::EchoData(in, out);

	/* echo parameters */
	in >> fGroupNumber;
	out << "\n Source element group. . . . . . . . . . . . . . = " << fGroupNumber << '\n';
	
	/* correct numbering */
	fGroupNumber--;
}

/* initialize data */
void MFPenaltySphereT::Initialize(void)
{
	/* inherited */
	PenaltySphereT::Initialize();
	
	/* allocate workspace */
	fCoords.Dimension(fNumContactNodes, rCoords.MinorDim());
	fCurrCoords.Dimension(fNumContactNodes, rEqnos.MinorDim());
	
	/* collect coordinates */
	fCoords.RowCollect(fContactNodes, rCoords);
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
	/* get pointer */
	fElementGroup = fFEManager.ElementGroup(fGroupNumber);
	if (!fElementGroup)
	{
		cout << "\n MFPenaltySphereT::SetElementGroup: error retrieving pointer\n"
		     <<   "     for group " << fGroupNumber + 1 << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* check */
	if (fElementGroup->InterpolantDOFs())
	{
		cout << "\n MFPenaltySphereT::SetElementGroup: element group " << fGroupNumber + 1
		     << '\n' <<   "     has interpolant DOF's. Use PenaltySphereT." << endl;
		throw ExceptionT::kGeneralFail;
	}
}
