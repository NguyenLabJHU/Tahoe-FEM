/* $Id: FEManagerT_bridging.cpp,v 1.1.2.1 2003-02-06 02:39:08 paklein Exp $ */
#include "FEManagerT_bridging.h"
#include "ModelManagerT.h"

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm):
	FEManagerT(input, output, comm)
{

}

/* initialize nodes that follow the field computed by this instance */
void SetFollowers(const iArrayT& nodes, const StringT& field, NodeManagerT& node_manager)
{
	// map incoming nodes into elements
	
	// store local nodes and shape function weights to do the interpolation

}

/* initialize data for the driving field */
void SetExactSolution(const iArrayT& nodes, const StringT& field,  NodeManagerT& node_manager)
{
	// map incoming nodes into elements
	
	// mark affected nodes as prescribed
	
	// reset the group equations numbers
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* initialize members */
void FEManagerT_bridging::ReadParameters(InitCodeT init)
{
	/* inherited */
	FEManagerT::ReadParameters(init);

	/* collect ghost nodes */
	ArrayT<StringT> id_list;
	fModelManager->NodeSetList(fMainIn, id_list);
	fModelManager->ManyNodeSets(id_list, fGhostNodes);
}

/* write parameters to main out */
void FEManagerT_bridging::WriteParameters(void) const
{
	/* inherited */
	FEManagerT::WriteParameters();

	/* echo ghost nodes */
	fMainOut << "\n Ghost nodes:\n";
	fMainOut << fGhostNodes.wrap(10) << '\n';
}

/* construct node manager */
void FEManagerT_bridging::SetNodeManager(void)
{
	/* inherited */
	FEManagerT::SetNodeManager();
	
	//declare ghost nodes with node manager
}
