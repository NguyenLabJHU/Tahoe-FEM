/* $Id: FEManagerT_bridging.cpp,v 1.1.2.2 2003-02-10 02:22:04 paklein Exp $ */
#include "FEManagerT_bridging.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "BridgingScaleT.h"
#include "KBC_PrescribedT.h"
#include "KBC_CardT.h"

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT(input, output, comm),
	fBridgingIn(bridging_input),
	fBridgingScale(NULL),
	fSolutionDriver(NULL)
{

}

/* initialize the ghost node information */
void FEManagerT_bridging::InitGhostNodes(void)
{
	const char caller[] = "FEManagerT_bridging::InitGhostNodes";

	/* collect ghost nodes */
	ArrayT<StringT> id_list;
	fModelManager->NodeSetList(fBridgingIn, id_list);
	fModelManager->ManyNodeSets(id_list, fGhostNodes);

	/* assume atomistic field is "displacement" */
	StringT field = "displacement";
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* create controller to driver solution */
	if (!fSolutionDriver) {
	
		/* construct new contoller */
		fSolutionDriver = new KBC_PrescribedT(*fNodeManager);
	
		/* add to field */
		the_field->AddKBCController(fSolutionDriver);
	}

	/* generate KBC cards - all degrees of freedom */
	int ndof = the_field->NumDOF();
	ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
	KBC_cards.Dimension(fGhostNodes.Length()*ndof);
	int dex = 0;
	for (int j = 0; j < ndof; j++)
		for (int i = 0; i < fDrivenCellNodes.Length(); i++)
			KBC_cards[dex++].SetValues(fGhostNodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	
	/* reset the group equations numbers */
	SetEquationSystem(the_field->Group());

	/* echo ghost nodes */
	if (fPrintInput) {
		fMainOut << "\n Ghost nodes:\n";
		fMainOut << fGhostNodes.wrap(10) << '\n';
	}

	/* mark nodes as ghost */
	fNonGhostNodes.Dimension(fModelManager->NumNodes() - fGhostNodes.Length());
	iArrayT is_ghost(fModelManager->NumNodes());
	is_ghost = 0;
	for (int i = 0; i < fGhostNodes.Length(); i++)
		is_ghost[fGhostNodes[i]] = 1;

	/* check for uniqueness */
	int ng = is_ghost.Count(1);
	if (ng != fGhostNodes.Length())
		ExceptionT::GeneralFail(caller, "list of ghost nodes contains %d duplicates",
			fGhostNodes.Length() - ng);

	/* collect non-ghost nodes */
	dex = 0;
	for (int i = 0; i < is_ghost.Length(); i++)
		if (is_ghost[i] == 0)
			fNonGhostNodes[dex++] = i;
}

/* initialize nodes that follow the field computed by this instance */
void FEManagerT_bridging::InitInterpolation(const iArrayT& nodes, const StringT& field, NodeManagerT& node_manager)
{
#pragma unused(field)

	/* map nodes into cells (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().MaptoCells(nodes, &init_coords, NULL, fFollowerCellData);
}

/* initialize data for the driving field */
void FEManagerT_bridging::InitProjection(const iArrayT& nodes, const StringT& field,  NodeManagerT& node_manager)
{
	const char caller[] = "FEManagerT_bridging::SetExactSolution";

	/* map nodes into cells (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().MaptoCells(nodes, &init_coords, NULL, fDrivenCellData);

	/* get the associated field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* create controller to driver solution */
	if (!fSolutionDriver) {
	
		/* construct new contoller */
		fSolutionDriver = new KBC_PrescribedT(*fNodeManager);
	
		/* add to field */
		the_field->AddKBCController(fSolutionDriver);
	}

	/* collect nodes in non-empty cells */
	fDrivenCellData.CollectCellNodes(fDrivenCellNodes);

	/* generate KBC cards - all degrees of freedom */
	int ndof = the_field->NumDOF();
	ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
	KBC_cards.Dimension(fDrivenCellNodes.Length()*ndof);
	int dex = 0;
	for (int j = 0; j < ndof; j++)
		for (int i = 0; i < fDrivenCellNodes.Length(); i++)
			KBC_cards[dex++].SetValues(fDrivenCellNodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	
	/* reset the group equations numbers */
	SetEquationSystem(the_field->Group());
}

/*************************************************************************
 * Private
 *************************************************************************/

/* the bridging scale element group */
BridgingScaleT& FEManagerT_bridging::BridgingScale(void) const
{
	/* find bridging scale group */
	if (!fBridgingScale) {
	
		/* search through element groups */
		for (int i = 0; !fBridgingScale && i < fElementGroups.Length(); i++)
		{
			/* try cast */
			ElementBaseT* element_base = fElementGroups[i];
			
			/* need non-const pointer to this */
			FEManagerT_bridging* fe = (FEManagerT_bridging*) this;
			fe->fBridgingScale = dynamic_cast<BridgingScaleT*>(element_base);
		}
		
		/* not found */
		if (!fBridgingScale)
			ExceptionT::GeneralFail("FEManagerT_bridging::BridgingScale",
				"did not find BridgingScaleT element group");
	}
	
	return *fBridgingScale;
}
