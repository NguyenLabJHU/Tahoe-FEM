/* $Id: SetOfNodesKBCT.cpp,v 1.4 2004-07-15 08:31:21 paklein Exp $ */
#include "SetOfNodesKBCT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* constructor */
SetOfNodesKBCT::SetOfNodesKBCT(const BasicSupportT& support, BasicFieldT& field):
	KBC_ControllerT(support),
	fSchedule(NULL)
{
#pragma unused(field)
	// nein
}

SetOfNodesKBCT::~SetOfNodesKBCT(void)
{
   
}

void SetOfNodesKBCT::Initialize(ifstreamT& in)
{
	int readRegion; // 0 for set, 1 for Region
	in >> readRegion;
	
	if (!readRegion)  // Initialize a range of nodes
	{		
		int iMin, iMax; // Range of Nodes to control
		
		in >> iMin >> iMax; 
		if (iMin > iMax)
			ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Bad node ID range specification \n");
		iMin--;
		iMax--;
	
		int numBCs;
		in >> numBCs;
		if (numBCs <= 0)
			ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Expecting numBCs > 0");
	
		iArrayT dof(numBCs), fScheduleNum(numBCs);
		dArrayT fScale(numBCs);
		ArrayT<KBC_CardT::CodeT> code(numBCs);
	 
	 	for (int i = 0; i < numBCs; i++)
	 	{
			in >> dof[i];
			dof[i]--;
	
			int i_code;
			in >> i_code;
			code[i] = KBC_CardT::int2CodeT(i_code);
	
			in >> fScheduleNum[i] >> fScale[i]; 
		
			if (code[i] != KBC_CardT::kFix)
			{
				fScheduleNum[i]--;
				fSchedule = fSupport.Schedule(fScheduleNum[i]);	
				if (!fSchedule) 
					ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Cannot get schedule %d \n",fScheduleNum[i]);
			}
			
		}
		
		fKBC_Cards.Dimension((iMax - iMin + 1) * numBCs);

		/* generate BC cards */
		KBC_CardT* pcard = fKBC_Cards.Pointer();
		
		for (int i = 0; i < numBCs; i++)
		{
			for (int j = 0; j < iMax - iMin + 1; j++, pcard++)
			{			
				/* set values */
				const ScheduleT* schedule = fSupport.Schedule(fScheduleNum[i]);
				pcard->SetValues(iMin + j, dof[i], code[i], schedule, fScale[i]);
			}
	
		}
		
	}
	else // Initialize regions
	{		
		int numSD = fSupport.NumSD();
		dArrayT fxmin(numSD), fxmax(numSD);
		
		/* get nodal coordinates */
		const dArray2DT& coords = fSupport.InitialCoordinates();

		in >> fxmin;
		in >> fxmax;
		in >> nIncs; 
		if (nIncs < 0) 
			ExceptionT::BadInputValue("SetOfNodesKBCT::InitRegion","Bad increment value");

		for (int j = 0; j < numSD; j++)
			if (fxmin[j] >= fxmax[j])
				ExceptionT::BadInputValue("SetOfNodesKBCT::InitRegion","Bad bounding box coordinates");
	
		if (fxmin.Length() != coords.MinorDim())
		ExceptionT::GeneralFail("ThermostattedRegionT::NodesInRegion",
				"Dimension mismatch between coords and bounding box");

		AutoArrayT<int> tmpList;		

		double* xmin = fxmin.Pointer();
		double* xmax = fxmax.Pointer(); 
		const double* x_j;
		int ihits = 0;
		int nnd = coords.MajorDim();
		for (int j = 0; j < nnd; j++)
		{
			bool inBox = true;
			x_j = coords(j);
			for (int k = 0; inBox && k < numSD; k++)
			{
				inBox = (xmin[k] < x_j[k]) && (x_j[k] < xmax[k]);
			}
			if (inBox)
			{
				tmpList.Append(j);
				ihits++;
			}
		}
			
		fNodes.Dimension(ihits);
		tmpList.CopyInto(fNodes);
		
		// save a bit of memory just in case
		tmpList.Free();
		
		int numBCs;
		in >> numBCs;
		if (numBCs <= 0)
			ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Expecting numBCs > 0");
	
		iArrayT dof(numBCs), fScheduleNum(numBCs);
		dArrayT fScale(numBCs);
		ArrayT<KBC_CardT::CodeT> code(numBCs);
	 
	 	for (int i = 0; i < numBCs; i++)
	 	{
			in >> dof[i];
			dof[i]--;
	
			int i_code;
			in >> i_code;
			code[i] = KBC_CardT::int2CodeT(i_code);
	
			in >> fScheduleNum[i] >> fScale[i]; 
		
			if (code[i] != KBC_CardT::kFix)
			{
				fScheduleNum[i]--;
				fSchedule = fSupport.Schedule(fScheduleNum[i]);	
				if (!fSchedule) 
					ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Cannot get schedule %d \n",fScheduleNum[i]);
			}
			
		}
		
		fKBC_Cards.Dimension(fNodes.Length()*numBCs);

		/* generate BC cards */
		KBC_CardT* pcard = fKBC_Cards.Pointer();
		
		for (int i = 0; i < numBCs; i++)
		{
			// loop over number of KBCs per node set
		
			for (int j = 0; j < fNodes.Length(); j++, pcard++)
			{			
				/* set values */
				const ScheduleT* schedule = fSupport.Schedule(fScheduleNum[i]);
				pcard->SetValues(fNodes[j], dof[i], code[i], schedule, fScale[i]);
			}
		
		}
	}	
}



