/* $Id: PenaltyWallT.cpp,v 1.9 2003-01-29 07:35:22 paklein Exp $ */
/* created: paklein (02/25/1997)                                          */

#include "PenaltyWallT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "fstreamT.h"
#include "FEManagerT.h"
#include "eIntegratorT.h"
#include "Vector3T.h"


using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
PenaltyWallT::PenaltyWallT(FEManagerT& fe_manager,
	int group,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	PenaltyRegionT(fe_manager, group, eqnos, coords, vels),

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
void PenaltyWallT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{	
#pragma unused(sys_type)
#if 0
	//TEMP
	if (fmu > kSmall)
	{
		cout << "\n PenaltyWallT::ApplyLHS: tangent is frictionless only" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
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
			fFEManager.AssembleLHS(fGroup, fLHS, fi_sh);
		}
	}
}

/* input processing */
void PenaltyWallT::EchoData(ifstreamT& in, ostream& out)
{
	/* inherited */
	PenaltyRegionT::EchoData(in, out);

	/* echo parameters */
	in >> fnormal; fnormal.UnitVector();
	out << " Wall normal:\n" << fnormal << '\n';
	
	/* transformation tensor */
	if (rCoords.MinorDim() == 2)
	{
		/* set column vectors */
		fQ.SetCol(0, fnormal);
		fQ(0,1) =-fnormal[1];
		fQ(1,1) = fnormal[0];
	}
	else if (rCoords.MinorDim() == 3)
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
	else throw ExceptionT::kGeneralFail;
}

/* initialize data */
void PenaltyWallT::Initialize(void)
{
	/* inherited */
	PenaltyRegionT::Initialize();

	/* memory relative displacements */
	fp_i.Dimension(fNumContactNodes, rCoords.MinorDim());
	fv_i.Dimension(fNumContactNodes, rCoords.MinorDim());	
}

/**********************************************************************
* Private
**********************************************************************/

/* compute the nodal contribution to the residual force vector */
void PenaltyWallT::ComputeContactForce(double kforce)
{
	/* with "friction */
//	if (fmu > kSmall)
	if (false)
	{
		//TEMP
		cout << "\n PenaltyWallT::ComputeContactForce: general (2D/3D) friction implementation\n";
		cout <<   "     is not available" << endl;
		throw ExceptionT::kGeneralFail;

		if (!pVels)
		{
			cout << "\n PenaltyWallT::ComputeContactForce: velocities required with friction";
			cout << endl;
			throw ExceptionT::kGeneralFail;
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
