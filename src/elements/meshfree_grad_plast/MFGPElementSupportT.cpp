/* $Id: MFGPElementSupportT.cpp,v 1.1 2005-04-26 22:32:28 kyonten Exp $ */
#include "MFGPElementSupportT.h"

#include "iAutoArrayT.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "ElementCardT.h"
#include "D3MeshFreeSupportT.h"
#include "ElementBaseT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"

#include "ModelManagerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* parameters */
const int kHeadRoom = 10; // percent

/* constructor */
MFGPElementSupportT::MFGPElementSupportT(void):
	//ParameterInterfaceT("mfgp_element_support"),
	fShapes_displ(NULL),
	fShapes_plast(NULL),
	fElemNodesEXDispl(NULL),
	fElemNodesEXPlast(NULL),
	fNumElemNodesDispl(0),
	fNumElemNodesPlast(0),
	fNEEMatrix(kHeadRoom, true)
{

}

/* accessors */
D3MeshFreeSupportT& MFGPElementSupportT::D3MeshFreeSupport(void) const {
	return fShapes_displ->D3MeshFreeSupport();
}

void MFGPElementSupportT::Register(dMatrixT& m, LocalArrayT& field_1, LocalArrayT& field_2 ) { 

	fLocField_1 = field_1;
	fLocField_2 = field_2;
	
	fNEEMatrix.Register(m);
}

/* information about subordinate parameter lists */
//void MFGPElementSupportT::DefineSubs(SubListT& sub_list) const
//{
	/* inherited */
	//ParameterInterfaceT::DefineSubs(sub_list);
//}

/* accept parameter list */
//void MFGPElementSupportT::TakeParameterList(const ParameterListT& list)
//{
	/* inherited */
	//ParameterInterfaceT::TakeParameterList(list);
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization */
void MFGPElementSupportT::InitSupport(ostream& out,
	AutoArrayT<ElementCardT>& elem_cards_displ, AutoArrayT<ElementCardT>& elem_cards_plast,
	const iArrayT& surface_nodes, int numDOF_displ, int numDOF_plast, 
	int max_node_num, ModelManagerT* model)
{
	/* configure variable length element arrays */
	fElemNodesEXDispl = &(fShapes_displ->ElementNeighbors());
	fElemEqnosEXDispl.Configure(fShapes_displ->ElementNeighborsCounts(), numDOF_displ);
	
	fElemNodesEXPlast = &(fShapes_plast->ElementNeighbors());
	fElemEqnosEXPlast.Configure(fShapes_plast->ElementNeighborsCounts(), numDOF_plast);
	
	/* set element card pointers */
	int num_cells = elem_cards_displ.Length();
	fUNodeLists.Dimension(num_cells);
	for (int i = 0; i < num_cells; i++)
	{
		ElementCardT& card = elem_cards_displ[i];
	
		/* field nodes */
		fElemNodesEXDispl->RowAlias(i, fUNodeLists[i]);
		card.SetNodesU(fUNodeLists[i]);

		/* field equations */
		fElemEqnosEXDispl.RowAlias(i, card.Equations());
	}
}

/* redimension off-diagonal stiffness matrices */
void MFGPElementSupportT::SetOffDiagMatrix(int element)
{
	/* current number of element neighbors */
	fNumElemNodesDispl = fElemNodesEXDispl->MinorDim(element);
	int neq_displ = fNumElemNodesDispl*fLocField_1.MinorDim();
	
	fNumElemNodesPlast = fElemNodesEXPlast->MinorDim(element);
	int neq_plast = fNumElemNodesPlast*fLocField_2.MinorDim();

	/* redimension workspace arrays */
	if (fLocField_1.MinorDim() > fLocField_2.MinorDim())
		fNEEMatrix.Dimension(neq_displ, neq_plast);
	else
		fNEEMatrix.Dimension(neq_plast, neq_displ);
}
