/* $Id: SymmetricNodesT.cpp,v 1.4.4.1 2002-12-10 17:08:52 paklein Exp $ */
#include "SymmetricNodesT.h"
#include "AutoArrayT.h"
#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"
#include "BasicFieldT.h"
#include "ElementsConfig.h"

#ifdef COHESIVE_SURFACE_ELEMENT
#include "TiedPotentialT.h"
#endif

using namespace Tahoe;

/* constructor */
SymmetricNodesT::SymmetricNodesT(NodeManagerT& node_manager, BasicFieldT& field):
	TiedNodesT(node_manager, field)
{
#ifndef COHESIVE_SURFACE_ELEMENT
	ExceptionT::BadInputValue("TiedNodesT::TiedNodesT", "COHESIVE_SURFACE_ELEMENT not enabled");
#endif
}

/* initialize data. Must be called immediately after construction */
void SymmetricNodesT::Initialize(ifstreamT& in)
{
	/* inherited */
	TiedNodesT::Initialize(in);

	int nsd = fNodeManager.NumSD();
	dArrayT fv1(nsd);	

	/* Read in points in symmetry plane */
	fDir.Dimension(nsd+1,nsd);
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
//void SymmetricNodesT::CopyKinematics(void)
//{
//	int nsd = fNodeManager.NumSD();
//	dArrayT w(nsd);
//	for (int i = 0; i < fPairStatus.Length(); i++)
//		if (fPairStatus[i] == kTied)
//		{
			/* destination and source */
//			int follower = fNodePairs(i,0);
//			int leader   = fNodePairs(i,1);

			/* kinematics */
//			for (int j = 0; j < fField.Order(); j++)
//			{
//				dArray2DT& u = fField[j];

				/* copy data from the leader */	
//				u.RowCopy(leader,w);

				/* modify data for reflection symmetry */
//				double scalarProd = 0.;
//				for (int k = 0; k < nsd; k++)
//				  scalarProd += w[k]*fDir(0,k);
//				scalarProd *= 2.;
//				for (int k = 0;k < nsd;k++)
//				  w[k] -= scalarProd*fDir(0,k);
//				u.SetRow(follower,w);
//			}
//		}
//}
bool SymmetricNodesT::ChangeStatus(void)
{
#ifdef COHESIVE_SURFACE_ELEMENT
  	bool changeQ = false;
	ElementBaseT* surroundingGroup = fFEManager.ElementGroup(0);
  	if (!surroundingGroup)
    {
      	cout <<" Group 0 doesn't exist \n";
      	throw ExceptionT::kGeneralFail;
    }
  	surroundingGroup->SendOutput(TiedPotentialT::kAverageCode);
  	dArray2DT fNodalQs = fNodeManager.OutputAverage();

  	for (int i = 0; i < fNodePairs.MajorDim();i++) 
    {		
  		dArrayT sigma(fNodalQs.MinorDim(),fNodalQs(fNodePairs(i,0)));
	  	if (fPairStatus[i] == kTied && TiedPotentialT::InitiationQ(sigma.Pointer()))     
		{ 
	  		fPairStatus[i] = kFree;
	  		changeQ = true;
	  		cout << "Freed a node !!! \n";
		}
    }

  	return changeQ;
#else
	return false;
#endif
}

