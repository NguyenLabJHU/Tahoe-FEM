/* $Id: AtomicConveyorT.cpp,v 1.1.2.2 2003-10-02 19:31:28 cjkimme Exp $ */
#include "AtomicConveyorT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "RandomNumberT.h"
#include "MessageT.h"
#include "CommunicatorT.h"
#include "iArrayT.h"
#include "ArrayT.h"
#include "StringT.h"
#include "KBC_CardT.h"
#include "ConveyorParticleT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
AtomicConveyorT::AtomicConveyorT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fField(field),
	fTopNodes(),
	fBottomNodes(),
	fSD(fNodeManager.NumSD())
{
	if (fSD != 2) ExceptionT::GeneralFail("AtomicConeyorT::AtomicConveyorT","Whoa hoss. This needs to be 2D for now\n");
}

AtomicConveyorT::~AtomicConveyorT(void)
{

}

void AtomicConveyorT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

}

void AtomicConveyorT::Initialize(ifstreamT& in)
{
	fBoxMin.Dimension(fSD);
	fBoxMax.Dimension(fSD);
	fBoxSize.Dimension(fSD);

	/* Get info about the crystal's reference (unstrained) configuration */
	in >> fBoxMin;
	in >> fBoxMax;
	
	for (int i = 0; i < fSD; i++)
		if (fBoxMin[i] > fBoxMax[i])
			ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bad box dimensions\n");
	
	for (int i = 0; i < fSD; i++) // useful quantity for making a new, unstrained crystal
		fBoxSize[i] = fBoxMax[i] - fBoxMin[i];
	
//	in >> initStrain; // This is an IC. Should it be moved to be with the rest of the ICs

	/* figure out how top/bottom BCs are specified */
	int BCformat; 
	in >> BCformat;
	
	if (!BCformat) // 0 means node set IDs will be given
	{
		ArrayT<StringT> fNodeIds;
		
		ReadNodes(in, fNodeIds, fTopNodes);
		
		if (fTopNodes.Length() == 0)
			ExceptionT::GeneralFail("AtomicConveyorT::Initialize","No top node sets to apply BCs to\n");
			
		ReadNodes(in, fNodeIds, fBottomNodes);
		
		if (fBottomNodes.Length() == 0)
			ExceptionT::GeneralFail("AtomicConveyorT::Initialize","No top node sets to apply BCs to\n");
		
	}
	else
	{
		fTopMin.Dimension(fSD);
		fTopMax.Dimension(fSD);
		fBottomMin.Dimension(fSD);
		fBottomMax.Dimension(fSD);
	
		if (BCformat == 1) // regions are given 
		{
			in >> fTopMin;
			in >> fTopMax;
			
			for (int i = 0; i < fSD; i++)
				if (fTopMin[i] >= fTopMax[i])
					ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bad top bounding box coordinates");		
			
			in >> fBottomMin;
			in >> fBottomMax;
			
			for (int i = 0; i < fSD; i++)
				if (fBottomMin[i] >= fBottomMax[i])
					ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bad bottom bounding box coordinates");		
					
			if (fTopMin[fSD-1] <= fBottomMax[fSD-1])
					ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bottom region is higher than top one");		
			
		}
		else // find nodes within tolerance of top and bottom of system 
		{
			double TopTol, BottomTol;
			in >> TopTol;
			in >> BottomTol;
			
			if (TopTol <= 0)
				ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Zero or negative width of top BC region\n");
			
			if (BottomTol <= 0)
				ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Zero or negative width of bottom BC region\n");
			
			fTopMin = fBoxMin;
			fTopMax = fBoxMax;
			fBottomMin = fBoxMin;
			fBottomMax = fBoxMax;
				
			fTopMin[fSD - 1] = fTopMax[fSD - 1] - TopTol;
			fBottomMax[fSD - 1] = fBottomMin[fSD - 1] + BottomTol;

		}

		NodesInRegion(fTopMin,fTopMax,fTopNodes);
		NodesInRegion(fBottomMin,fBottomMax,fBottomNodes);

	}
		
	/* get schedule for top and bottom BCs */
	int scheduleNum;
	in >> scheduleNum; scheduleNum--;
	in >> fScale; /* value for top of strip. value for bottom is negative of this */
	fSchedule = fNodeManager.Schedule(scheduleNum);	
	if (!fSchedule) throw ExceptionT::kBadInputValue;
	
	/* constrain the nodes just located with the values just read in */
	SetBCCards();
		
	/* read in left/right BCs. Ramped damping, cutting and pasting of new material */
	// width of left and right regions
	// maximum value of damping coefficients in each region
	
	// width of cut-and-paste region
	// how close crack tip gets to right edge before doing it	
		
	/* determine ICs */
	// need flag for minimization, lin elast solution, or what have you
	// perhaps also read coords from a file
	
	// initial velocities 
	// only near the crack tip or uniform strain rate
	
	/* crack tip specifications */
	in >> fPreCrackLength; 
	in >> fCrackHeight;
	in >> fTipRegionHeight;
	
	if ((fTipRegionHeight <= 0.) || (fTipRegionHeight > fBoxSize[fSD - 1]))
		ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bad Tip Region Height\n");
		
	int isDynamic;
	in >> isDynamic; 
	
	isDynamic ? qDynamic = true : qDynamic = false;
		
	// tolerance in locating precrack.
	
	in >> nTipIncs; // how often to locate the precrack
	if (nTipIncs < 0) ExceptionT::BadInputValue("AtomicConveyorT::Initialize",
				"Bad increment for locating crack tip given\n");
	
	in >> iElementGroup; 
	iElementGroup--;
	
	if (iElementGroup < 0)
		ExceptionT::BadInputValue("AtomicConveyorT::Initialize","Bad Element Block ID read\n");	

	
}

void AtomicConveyorT::InitStep(void)
{
					
	/* inherited */
	KBC_ControllerT::InitStep();

}

#pragma message("All my min/max for regions are based on 2D so far\n")
void AtomicConveyorT::InitialCondition(void)
{
	/* strain the system to the initial state */
	
	/* create precrack */
	dArrayT aboveRegionMin(fSD), aboveRegionMax(fSD), belowRegionMin(fSD), belowRegionMax(fSD);
	iArrayT NodesAboveCrack, NodesBelowCrack;
	
	aboveRegionMin = fBoxMin;
	aboveRegionMin[fSD - 1] = fCrackHeight;
	aboveRegionMax = fBoxMin; 
	aboveRegionMax[fSD - 2] += fPreCrackLength;
	aboveRegionMax[fSD - 1] = fCrackHeight + fTipRegionHeight;
	belowRegionMin = fBoxMin;
	belowRegionMin[fSD - 1] = fCrackHeight - fTipRegionHeight;
	belowRegionMax = fBoxMin;
	belowRegionMax[fSD - 2] += fPreCrackLength;
	belowRegionMax[fSD - 1] = fCrackHeight;
	
	NodesInRegion(aboveRegionMin, aboveRegionMax, NodesAboveCrack);
	NodesInRegion(belowRegionMin, belowRegionMax, NodesBelowCrack);
	
	cout << "AtomicConveyorT::InitialCondition "<< NodesAboveCrack.Length() << " above ";
	cout << NodesBelowCrack.Length() << " below \n";
	// turn off interactions between these nodes 
	
	ConveyorParticleT* elementGroup = dynamic_cast<ConveyorParticleT*>(fNodeManager.FEManager().ElementGroup(iElementGroup));
	
	if (!elementGroup)
		ExceptionT::GeneralFail("AtomicConveyorT::InitialCondition","Element group %d is not of type ConveyorParticleT",iElementGroup);

	elementGroup->CreateNoninteractingAtoms(NodesAboveCrack, NodesBelowCrack); 
	
	/* apply initial velocities */
	
	/* number of scaled nodes */
//	int n_scaled = 0;
//	int ndof = fField.NumDOF();
	
	/* get MPI stuff */
//	CommunicatorT& communicator = fNodeManager.FEManager().Communicator();
//	int nProcs = fNodeManager.Size();
//	int thisProc = fNodeManager.Rank();
//	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();

}

/**********************************************************************
 * Protected
 **********************************************************************/

void AtomicConveyorT::NodesInRegion(dArrayT& xmin, dArrayT& xmax, iArrayT& nodes)
{

	/* get nodal coordinates */
	const dArray2DT& coords = fNodeManager.InitialCoordinates();

	AutoArrayT<int> tmpList;
	double* x_i;
	int ihits = 0;
	int nnd = coords.MajorDim();
	for (int i = 0; i < nnd; i++)
	{
		bool inBox = true;
		x_i = coords(i);
		for (int j = 0; inBox && j < fSD; j++)
		{
			inBox = (xmin[j] <= x_i[j]) && (x_i[j] <= xmax[j]);
 		}
		if (inBox)
		{
			tmpList.Append(i);
			ihits++;
		}
	}
	nodes.Dimension(ihits);
	tmpList.CopyInto(nodes);
	tmpList.Free();
}

/* initialize the current step */
void AtomicConveyorT::SetBCCards(void)
{
	/* number of scaled nodes */
	int n_top = 0, n_bottom = 0; 
		
	/* get MPI stuff */
	CommunicatorT& communicator = fNodeManager.FEManager().Communicator();
	int nProcs = fNodeManager.Size();
	int thisProc = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();	

	/* figure out which nodes to affect */
	iArrayT myTopNodes, myBottomNodes;
	myTopNodes.Dimension(fTopNodes.Length());
	myBottomNodes.Dimension(fBottomNodes.Length());

	if (fNodeManager.Size() == 1 || !pMap)
	{
		myTopNodes.Set(fTopNodes.Length(), fTopNodes.Pointer());
		n_top = myTopNodes.Length();
		myBottomNodes.Set(fBottomNodes.Length(), fBottomNodes.Pointer());
		n_bottom = myBottomNodes.Length();
	}
	else
	{
		for (int i = 0; i < myTopNodes.Length(); i++)
			if ((*pMap)[fTopNodes[i]] == thisProc)
			{
				myTopNodes[i] = fTopNodes[i];
				n_top++;
			}
			else
				myTopNodes[i] = -1;
			
		for (int i = 0; i < myBottomNodes.Length(); i++)
			if ((*pMap)[fBottomNodes[i]] == thisProc)
			{
				myBottomNodes[i] = fBottomNodes[i];
				n_bottom++;
			}
			else
				myBottomNodes[i] = -1;
	}	
		
	/* generate BC cards */
	fKBC_Cards.Dimension(fSD*(n_top + n_bottom));
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < myTopNodes.Length(); i++)
	{
		if (myTopNodes[i] >= 0)
		{
			/* loading direction */
			pcard->SetValues(myTopNodes[i], fSD - 1, KBC_CardT::kDsp, 0, fScale);
			pcard->SetSchedule(fSchedule);
			pcard++;
			
			/* Fixed DOFs */
			pcard->SetValues(myTopNodes[i], fSD - 2, KBC_CardT::kFix, 0, 0.);
			pcard++;
			if (fSD == 3)
			{
				pcard->SetValues(myTopNodes[i], 0, KBC_CardT::kFix, 0, 0.);
				pcard++;
			}
			
		} 
	}
	for (int i = 0; i < myBottomNodes.Length(); i++)
	{
		if (myBottomNodes[i] >= 0)
		{
			/* loading direction */
			pcard->SetValues(myBottomNodes[i], fSD - 1, KBC_CardT::kDsp, 0, -fScale);
			pcard->SetSchedule(fSchedule);
			pcard++;
			
			/* Fixed DOFs */
			pcard->SetValues(myBottomNodes[i], fSD - 2, KBC_CardT::kFix, 0, 0.);
			pcard++;
			if (fSD == 3)
			{
				pcard->SetValues(myBottomNodes[i], 0, KBC_CardT::kFix, 0, 0.);
				pcard++;
			}
		} 
	}

}

