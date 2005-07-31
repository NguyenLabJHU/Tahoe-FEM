/* $Id: PenaltyWallT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (02/25/1997)                                          */

#include "PenaltyWallT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "fstreamT.h"
#include "FEManagerT.h"
#include "eControllerT.h"

const double Pi = acos(-1.0);

/* constructor */
PenaltyWallT::PenaltyWallT(FEManagerT& fe_manager,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	PenaltyRegionT(fe_manager, eqnos, coords, vels),

	/* wall normal and tangents */	
	fnormal(rCoords.MinorDim()),
	fntforce(rCoords.MinorDim()),
	fxyforce(rCoords.MinorDim()),
	fQ(rCoords.MinorDim()),
	
	/* work space */
	fLHS(eqnos.MinorDim(), ElementMatrixT::kSymmetric)
{

}

/* tangent */
void PenaltyWallT::ApplyLHS(void)
{	
	//TEMP
	if (fmu > kSmall)
	{
		cout << "\n PenaltyWallT::ApplyLHS: tangent is frictionless only" << endl;
		throw eGeneralFail;
	}

	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

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
			rEqnos.RowAlias(fContactNodes[i], fi_sh);
			fFEManager.AssembleLHS(fLHS, fi_sh);
		}
	}
}

/* input processing */
void PenaltyWallT::EchoData(ifstreamT& in, ostream& out)
{
	/* inherited */
	PenaltyRegionT::EchoData(in, out);

	/* echo parameters */
	in >> ftheta;
	in >> fmu;    if (fmu < 0.0) throw eBadInputValue;

	out << " Orientation of normal wrt x-axis (degrees). . . = " << ftheta << '\n';
	out << " Penalty stiffness . . . . . . . . . . . . . . . = " << fk << '\n';

	/* compute normal, tangent, and Q */
	ftheta *= Pi/180.0;
	
	if (rCoords.MinorDim() == 2)
	{
		fnormal[0] = cos(ftheta);
		fnormal[1] = sin(ftheta);

		fQ(0,0) = fQ(1,1) = cos(ftheta);
		fQ(1,0) = sin(ftheta);
		fQ(0,1) =-sin(ftheta);
	}
	else // 3D still only has rotation about z-axis
	{
		fnormal[0] = cos(ftheta);
		fnormal[1] = sin(ftheta);
		fnormal[2] = 0.0;

		fQ = 0.0;
		fQ(0,0) = fQ(1,1) = cos(ftheta);
		fQ(1,0) = sin(ftheta);
		fQ(0,1) =-sin(ftheta);
		fQ(2,2) = 1.0;
	}
}

/* initialize data */
void PenaltyWallT::Initialize(void)
{
	/* inherited */
	PenaltyRegionT::Initialize();

	/* memory relative displacements */
	fp_i.Allocate(fNumContactNodes, rCoords.MinorDim());
	fv_i.Allocate(fNumContactNodes, rCoords.MinorDim());	
}

/**********************************************************************
* Private
**********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltyWallT::ComputeContactForce(double kforce)
{
	/* with "friction */
	if (fmu > kSmall)
	{
		//TEMP
		cout << "\n PenaltyWallT::ComputeContactForce: general (2D/3D) friction implementation\n";
		cout <<   "     is not available" << endl;
		throw eGeneralFail;

		if (!pVels)
		{
			cout << "\n PenaltyWallT::ComputeContactForce: velocities required with friction";
			cout << endl;
			throw eGeneralFail;
		}
	
		/* compute relative positions and velocities */
		fp_i.RowCollect(fContactNodes, rCoords);
		fv_i.RowCollect(fContactNodes,*pVels);
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
				// old 2D implementation:
				// double tangent_comp = fv_i.DotRow(i, ftangent);	
				fntforce[0] =-fk*normal_comp*kforce;		
				//fntforce[1] = ((tangent_comp > 0.0) ? -1.0 : 1.0)*fmu*fntforce[0]*kforce;
		
				/* transform to x-y coordinates */
				fQ.Multx(fntforce, fxyforce);
				
				fContactForce2D.SetRow(i,fxyforce);
			}
		}
	}
	else
	{
		/* compute relative positions */
		fp_i.RowCollect(fContactNodes, rCoords);
		for (int j = 0; j < fNumContactNodes; j++)
			fp_i.AddToRowScaled(j, -1.0, fx);
	
		/* compute contact forces */
		fh_max = 0.0;
		fntforce = 0.0;
		fContactForce2D = 0.0;	
		for (int i = 0; i < fNumContactNodes; i++)
		{
			double normal_comp = fp_i.DotRow(i, fnormal);
		
			/* penetration */
			if (normal_comp < 0.0)
			{
				/* store max penetration */
				fh_max = (-normal_comp > fh_max) ? -normal_comp : fh_max;
			
				/* normal force */
				fntforce[0] =-fk*normal_comp*kforce;		
		
				/* transform to x-y coordinates */
				fQ.Multx(fntforce, fxyforce);
				
				fContactForce2D.SetRow(i, fxyforce);
			}
		}
	}
}
