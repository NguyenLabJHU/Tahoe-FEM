/* $Id: SymmetricNodesT.cpp,v 1.1.2.2 2002-04-30 01:30:21 paklein Exp $ */
#include "SymmetricNodesT.h"
#include "AutoArrayT.h"
#include "NodeManagerT.h"
//#include "ElementBaseT.h"
//#include "FEManagerT.h"
#include "ifstreamT.h"
#include "BasicFieldT.h"

/* constructor */
SymmetricNodesT::SymmetricNodesT(NodeManagerT& node_manager, BasicFieldT& field):
	/*KBC_ControllerT(node_manager),
	fEqnos(NULL),
	fKinematics(0),
	fDummySchedule(1.0),
	fFEManager(node_manager.FEManager())*/
	TiedNodesT(node_manager, field)
{

}

/* initialize data. Must be called immediately after construction */
void SymmetricNodesT::Initialize(ifstreamT& in)
{
	/* inherited */
	TiedNodesT::Initialize(in);

	int nsd = fNodeManager.NumSD();
	dArrayT fv1(nsd);	

	/* Read in points in symmetry plane */
	fDir.Allocate(nsd+1,nsd);
	for (int i = nsd;i > 1;i--)
	{
	    in >> fv1;
	    fDir.SetRow(i-1,fv1);
	}
	in >> fv1;	

	/* Last row contains the offset point */
	fDir.SetRow(nsd,fv1);

	/* second row through next-to-last one contain direction vectors */
	fv1 *= -1.;
	for (int i = 1;i < nsd; i++)
	{
	    fDir.AddToRowScaled(i,1.,fv1);
	}

	/* compute normal directions and store in first row */
	if (nsd == 2) 
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
//UNUSED
#if 0
	int curCoordsIndex = fKinematics.Length() - 1;
#endif

	int nsd = fNodeManager.NumSD();
	dArrayT w(nsd);
	for (int i = 0; i < fPairStatus.Length(); i++)
		if (fPairStatus[i] == kTied)
		{
			/* destination and source */
			int follower = fNodePairs(i,0);
			int leader   = fNodePairs(i,1);

			/* kinematics */
			for (int j = 0; j < fField.Order(); j++)
			{
				dArray2DT& u = fField[j];

				/* copy data from the leader */	
				u.RowCopy(leader,w);

				/*if (j != curCoordsIndex)
				  for (k = 0; k < fSD; k++)
				  w[k] -=  fDir(fSD,k);*/
				double scalarProd = 0.;
				for (int k = 0; k < nsd; k++)
				  scalarProd += w[k]*fDir(0,k);
				scalarProd *= 2.;
				for (int k = 0;k < nsd;k++)
				  w[k] -= scalarProd*fDir(0,k);
				/*if (j != curCoordsIndex)
				  for (k = 0; k < fSD; k++) 
				  w[k] += fDir(fSD,k);*/
				u.SetRow(follower,w);
			}
		}
}
