/* $Id: SetOfNodesKBCT.cpp,v 1.1 2003-05-28 17:44:52 cjkimme Exp $ */
#include "SetOfNodesKBCT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
SetOfNodesKBCT::SetOfNodesKBCT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fSchedule(NULL)
{
#pragma unused(field)
	// nein
}

SetOfNodesKBCT::~SetOfNodesKBCT(void)
{
   
}

void SetOfNodesKBCT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

}

void SetOfNodesKBCT::Initialize(ifstreamT& in)
{
	int numBCs;
	in >> numBCs;
	
	/* Read in data first */
	iArrayT iMin(numBCs), iMax(numBCs), dof(numBCs), fScheduleNum(numBCs);
	dArrayT fScale(numBCs);
	ArrayT<KBC_CardT::CodeT> code(numBCs);
	int numCards = 0;
	for (int i = 0; i < numBCs; i++)
	{
//	int SetOrRegion; 
		in >> iMin[i] >> iMax[i]; 
		if (iMin[i] > iMax[i])
			ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Bad node ID range specification for set %d\n",i);
		iMin[i]--;
		iMax[i]--;
	 
		in >> dof[i];
		dof[i]--;
	
		in >> code[i];
	
		in >> fScheduleNum[i] >> fScale[i]; 
		
		if (code[i] != KBC_CardT::kFix)
		{
			fScheduleNum[i]--;
			fSchedule = fNodeManager.Schedule(fScheduleNum[i]);	
			if (!fSchedule) 
				ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Cannot get schedule %d \n",fScheduleNum[i]);
		}
			
		numCards += iMax[i] - iMin[i] + 1;
	}
	
	fKBC_Cards.Dimension(numCards);

	/* generate BC cards */
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	
	for (int i = 0; i < numBCs; i++)
	{
		for (int j = 0; j < iMax[i] - iMin[i] + 1; j++, pcard++)
		{			
			/* set values */
			pcard->SetValues(iMin[i] + j, dof[i], code[i], fScheduleNum[i], fScale[i]);
		}
		
	}
	
}



