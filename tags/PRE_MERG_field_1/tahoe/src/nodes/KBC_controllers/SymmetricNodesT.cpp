/* $Id: SymmetricNodesT.cpp,v 1.1 2002-04-19 21:12:58 cjkimme Exp $ */
#include "SymmetricNodesT.h"
#include "AutoArrayT.h"
#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"


/* constructor */
SymmetricNodesT::SymmetricNodesT(NodeManagerT& node_manager):
	/*KBC_ControllerT(node_manager),
	fEqnos(NULL),
	fKinematics(0),
	fDummySchedule(1.0),
	fFEManager(node_manager.FEManager())*/
        TiedNodesT(node_manager),
	fSD(node_manager.NumSD())
{

}

/* initialize data. Must be called immediately after construction */
void SymmetricNodesT::Initialize(ifstreamT& in)
{
	/* inherited */
	TiedNodesT::Initialize(in);
		
	dArrayT fv1(fSD);	

	/* Read in points in symmetry plane */
	fDir.Allocate(fSD+1,fSD);
	for (int i = fSD;i > 1;i--)
	{
	    in >> fv1;
	    fDir.SetRow(i-1,fv1);
	}
	in >> fv1;	

	/* Last row contains the offset point */
	fDir.SetRow(fSD,fv1);

	/* second row through next-to-last one contain direction vectors */
	fv1 *= -1.;
	for (int i = 1;i < fSD; i++)
	{
	    fDir.AddToRowScaled(i,1.,fv1);
	}

	/* compute normal directions and store in first row */
	if (fSD == 2) 
	{
	    fDir(0,0) = -fDir(1,1);
	    fDir(0,1) = fDir(1,0);
	    fDir.ScaleRow(0,1./sqrt(fDir(0,0)*fDir(0,0)+fDir(0,1)*fDir(0,1)));
	}
	else
	{
	    double *v = fDir(1);
	    double *w = fDir(2);
	    fDir(0,0) = v[1]*w[2]-v[2]*w[1];
	    fDir(0,1) = v[2]*w[0]-v[0]*w[2];
	    fDir(0,2) = v[0]*w[1]-v[1]*w[0];
	    fDir.ScaleRow(0,1./sqrt(fDir(0,0)*fDir(0,0)+fDir(0,1)*fDir(0,1)+fDir(0,2)*fDir(0,2)));
	}
}	

/* copy kinematic information from the leader nodes to the follower nodes */
void SymmetricNodesT::CopyKinematics(void)
{

        int curCoordsIndex = fKinematics.Length() - 1;
	double scalarProd;
	dArrayT w(fSD);
	int k;
	
	for (int i = 0; i < fPairStatus.Length(); i++)
		if (fPairStatus[i] == kTied)
		{
			/* destination and source */
			int follower = fNodePairs(i,0);
			int leader   = fNodePairs(i,1);

			/* kinematics */
			for (int j = 0; j < fKinematics.Length(); j++)
			{
				dArray2DT& u = *(fKinematics[j]);

				/* copy data from the leader */	
				u.RowCopy(leader,w);

				/*if (j != curCoordsIndex)
				  for (k = 0; k < fSD; k++)
				  w[k] -=  fDir(fSD,k);*/
				scalarProd = 0.;
				for (k = 0; k < fSD; k++)
				  scalarProd += w[k]*fDir(0,k);
				scalarProd *= 2.;
				for (k = 0;k < fSD;k++)
				  w[k] -= scalarProd*fDir(0,k);
				/*if (j != curCoordsIndex)
				  for (k = 0; k < fSD; k++) 
				  w[k] += fDir(fSD,k);*/
				u.SetRow(follower,w);
			}
		}
}
