/* $Id: TotalLagrangianCBSurfaceT.cpp,v 1.1 2005-06-29 17:39:46 paklein Exp $ */
#include "TotalLagrangianCBSurfaceT.h"

#include "ModelManagerT.h"

using namespace Tahoe;

/* constructor */
TotalLagrangianCBSurfaceT::TotalLagrangianCBSurfaceT(const ElementSupportT& support):
	TotalLagrangianT(support)
{
	SetName("total_lagrangian_CBsurface");
}

/* accept parameter list */
void TotalLagrangianCBSurfaceT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	TotalLagrangianT::TakeParameterList(list);
	
	/* collect surface element information */
	ArrayT<StringT> block_ID;
	ElementBlockIDs(block_ID);
	ModelManagerT& model_manager = ElementSupport().ModelManager();
	model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);
}

