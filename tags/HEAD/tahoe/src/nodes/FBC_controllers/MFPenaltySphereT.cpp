/* $Id: MFPenaltySphereT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/17/2000)                                          */

#include "MFPenaltySphereT.h"

#include <iostream.h>

#include "Constants.h"
#include "FEManagerT.h"
#include "ElementBaseT.h"
#include "fstreamT.h"

/* constructor */
MFPenaltySphereT::MFPenaltySphereT(FEManagerT& fe_manager,
	const iArray2DT& eqnos, const dArray2DT& coords, const dArray2DT* vels):
	PenaltySphereT(fe_manager, eqnos, coords, vels),
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
	fCoords.Allocate(fNumContactNodes, rCoords.MinorDim());
	fCurrCoords.Allocate(fNumContactNodes, rEqnos.MinorDim());
	
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
	fh_max = 0.0;
	fContactForce2D = 0.0;	
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* center to striker */
		fCurrCoords.RowCopy(i, fv_OP);
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
		throw eGeneralFail;
	}
	
	/* check */
	if (fElementGroup->InterpolantDOFs())
	{
		cout << "\n MFPenaltySphereT::SetElementGroup: element group " << fGroupNumber + 1
		     << '\n' <<   "     has interpolant DOF's. Use PenaltySphereT." << endl;
		throw eGeneralFail;
	}
}
