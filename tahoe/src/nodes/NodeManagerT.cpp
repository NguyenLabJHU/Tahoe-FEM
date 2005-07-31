/* $Id: NodeManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
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
int NodeManagerT::DegreesOfFreedom(int nsd) const { return nsd; }

