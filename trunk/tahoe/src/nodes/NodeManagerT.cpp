/* $Id: NodeManagerT.cpp,v 1.2 2001-06-29 23:51:41 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */
/* Field variables plus averging                                          */

#include "NodeManagerT.h"
#include "FEManagerT.h"

/* kinematic BC controllers */
#include "K_FieldT.h"
#include "BimaterialK_FieldT.h"
#include "MappedPeriodicT.h"

/* constructor */
NodeManagerT::NodeManagerT(FEManagerT& fe_manager):
	NodeManagerPrimitive(fe_manager)
{

}

/* allocate global nodal arrays */
void NodeManagerT::AllocateGlobal(void)
{
	/* inherited */
	NodeManagerPrimitive::AllocateGlobal();

	/* set nodal averaging space */
	SetNumAverageRows(fNumNodes);
}

/**********************************************************************
* Protected
**********************************************************************/

/* BC Controllers */
KBC_ControllerT* NodeManagerT::NewKBC_Controller(int code)
{
	switch(code)
	{
		case KBC_ControllerT::kK_Field:
			return new K_FieldT(*this);

		case KBC_ControllerT::kBimaterialK_Field:	
			return new BimaterialK_FieldT(*this);

		case KBC_ControllerT::kMappedPeriodic:	
			return new MappedPeriodicT(*this);

		default:
			/* inherited */
			return NodeManagerPrimitive::NewKBC_Controller(code);
	}
}

/* return the number of DOF's */
int NodeManagerT::DegreesOfFreedom(int nsd) const 
{ 
	/* select based on the analysis code */
	switch (fFEManager.Analysis())
	{
		case GlobalT::kLinStatic:
		case GlobalT::kLinDynamic:
		case GlobalT::kNLStatic:
		case GlobalT::kNLDynamic:
		case GlobalT::kDR:
		case GlobalT::kLinExpDynamic:
		case GlobalT::kNLExpDynamic:
		case GlobalT::kVarNodeNLStatic:
		case GlobalT::kVarNodeNLExpDyn:
		case GlobalT::kAugLagStatic:
			return nsd;
		
		case GlobalT::kLinStaticHeat:
		case GlobalT::kLinTransHeat:
			return 1;
			
		default:
			cout << "\n NodeManagerT::DegreesOfFreedom: could not resolve analysis code: " 
			     << fFEManager.Analysis() << endl;
			throw eGeneralFail;
	}	
	return 0; 	
}

