/* $Id: MFGPElementSupportT.cpp,v 1.3 2005-05-03 20:12:35 kyonten Exp $ */
#include "MFGPElementSupportT.h"


#include "iAutoArrayT.h"
#include "MeshFreeShapeFunctionT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "ElementCardT.h"
#include "MeshFreeSupportT.h"
#include "ElementBaseT.h"
#include "ParameterUtils.h"

#include "ModelManagerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* parameters */
const int kHeadRoom = 10; // percent

/* constructor */
MFGPElementSupportT::MFGPElementSupportT(void):
	ParameterInterfaceT("mfgp_element_support"),
	fMFShapes(NULL),
	fMFShapes_displ(NULL),
	fMFShapes_plast(NULL),
	fNodalShapes(NULL),
	fElemNodesEX(NULL),
	fElemNodesEX_displ(NULL),
	fElemNodesEX_plast(NULL),
	fLocGroup(kHeadRoom),
	fNumElemenNodes(0),
	fNumElemenNodes_displ(0),
	fNumElemenNodes_plast(0),
	fNEEArray(kHeadRoom, true),
	fNEEMatrix(kHeadRoom, true),
	fNEEMatrixOD(kHeadRoom, true),
	fFieldSet(false),
	fMapShift(-1),
	fOffGridID(NULL),
	fInterpolantID(NULL),
	fMeshlessID(NULL)
{

}

/* accessors */
MeshFreeSupportT& MFGPElementSupportT::MeshFreeSupport(void) const {
	return fMFShapes->MeshFreeSupport();
}

/* register non-square matrix */
void MFGPElementSupportT::Register(dMatrixT& m, LocalArrayT& field_1, LocalArrayT& field_2 ) { 

	fLocField_1 = field_1;
	fLocField_2 = field_2;
	
	fNEEMatrixOD.Register(m);
}

/* information about subordinate parameter lists */
void MFGPElementSupportT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* ID's for off grid nodes */
	sub_list.AddSub("off_grid_node_ID_list", ParameterListT::ZeroOrOnce);

	/* ID's for interpolant nodes */
	sub_list.AddSub("interpolant_node_ID_list", ParameterListT::ZeroOrOnce);

	/* ID's for forces meshless nodes */
	sub_list.AddSub("meshfree_node_ID_list", ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void MFGPElementSupportT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* collect pointers to class parameters */
	fOffGridID = list.List("off_grid_node_ID_list");
	fInterpolantID = list.List("interpolant_node_ID_list");
	fMeshlessID = list.List("meshfree_node_ID_list");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization */
void MFGPElementSupportT::InitSupport(ostream& out,
	AutoArrayT<ElementCardT>& elem_cards, const iArrayT& surface_nodes,
	int numDOF, int max_node_num, ModelManagerT* model)
{
	/* configure variable length element arrays */
	fElemNodesEX = &(fMFShapes->ElementNeighbors());
	fElemEqnosEX.Configure(fMFShapes->ElementNeighborsCounts(), numDOF);
	
	/* set element card pointers */
	int num_cells = elem_cards.Length();
	fUNodeLists.Dimension(num_cells);
	for (int i = 0; i < num_cells; i++)
	{
		ElementCardT& card = elem_cards[i];
	
		/* field nodes */
		fElemNodesEX->RowAlias(i, fUNodeLists[i]);
		card.SetNodesU(fUNodeLists[i]);

		/* field equations */
		fElemEqnosEX.RowAlias(i, card.Equations());
	}
	
	/* collect FE/meshfree nodes */
	CollectNodesData(out, max_node_num, model);

	/* collect interpolant nodes = (auto) + (FE) - (EFG) */
	SetAllFENodes(surface_nodes);

	out << " Final number of interpolant nodes . . . . . . . = ";
	out << fAllFENodes.Length() << '\n';
	if (fAllFENodes.Length() > 0)
	{
		/* correct offset for output */
		fAllFENodes++;
		out << fAllFENodes.wrap(6) << '\n';
		fAllFENodes--;
	}
}

/* initialization */
void MFGPElementSupportT::InitSupport(ostream& out,
	AutoArrayT<ElementCardT>& elem_cards_displ, AutoArrayT<ElementCardT>& elem_cards_plast, 
	const iArrayT& surface_nodes, int numDOF_displ, int numDOF_plast, int max_node_num, 
	ModelManagerT* model)
{
	/* configure variable length element arrays */
	fElemNodesEX_displ = &(fMFShapes_displ->ElementNeighbors());
	fElemEqnosEX_displ.Configure(fMFShapes_displ->ElementNeighborsCounts(), numDOF_displ);
	
	fElemNodesEX_plast = &(fMFShapes_plast->ElementNeighbors());
	fElemEqnosEX_plast.Configure(fMFShapes_plast->ElementNeighborsCounts(), numDOF_plast);
	
	/* set element card pointers */
	int num_cells = elem_cards_displ.Length();
	fUNodeLists_displ.Dimension(num_cells);
	fUNodeLists_plast.Dimension(num_cells);
	for (int i = 0; i < num_cells; i++)
	{
		ElementCardT& card_displ = elem_cards_displ[i];
		ElementCardT& card_plast = elem_cards_plast[i];
	
		/* field nodes */
		fElemNodesEX_displ->RowAlias(i, fUNodeLists_displ[i]);
		card_displ.SetNodesU(fUNodeLists_displ[i]);
		fElemNodesEX_plast->RowAlias(i, fUNodeLists_plast[i]);
		card_plast.SetNodesU(fUNodeLists_plast[i]);

		/* field equations */
		fElemEqnosEX_displ.RowAlias(i, card_displ.Equations());
		fElemEqnosEX_plast.RowAlias(i, card_plast.Equations());
		
		/* combine the two equations */
		//fElemEqnosEX_comb. ?? column copy
	}
	
	/* collect FE/meshfree nodes */
	CollectNodesData(out, max_node_num, model);

	/* collect interpolant nodes = (auto) + (FE) - (EFG) */
	SetAllFENodes(surface_nodes);

	out << " Final number of interpolant nodes . . . . . . . = ";
	out << fAllFENodes.Length() << '\n';
	if (fAllFENodes.Length() > 0)
	{
		/* correct offset for output */
		fAllFENodes++;
		out << fAllFENodes.wrap(6) << '\n';
		fAllFENodes--;
	}
}

/* resize number of field nodes */
int MFGPElementSupportT::SetElementNodes(int element)
{
	/* current number of element neighbors */
	fNumElemenNodes = fElemNodesEX->MinorDim(element);
	int neq = fNumElemenNodes*fLocGroup.MinorDim();
	
	/* redimension local arrays */
	fLocGroup.SetNumberOfNodes(fNumElemenNodes);

	/* redimension workspace arrays */
	fNEEArray.Dimension(neq, false);
	fNEEMatrix.Dimension(neq, neq);

	return fNumElemenNodes;
}

/* redimension off-diagonal stiffness matrices */
void MFGPElementSupportT::SetOffDiagMatrix(int element)
{
	/* current number of element neighbors */
	fNumElemenNodes_displ = fElemNodesEX_displ->MinorDim(element);
	int neq_displ = fNumElemenNodes_displ*fLocField_1.MinorDim();
	
	fNumElemenNodes_plast = fElemNodesEX_plast->MinorDim(element);
	int neq_plast = fNumElemenNodes_plast*fLocField_2.MinorDim();

	/* redimension workspace arrays */
	if (fLocField_1.MinorDim() > fLocField_2.MinorDim())
		fNEEMatrixOD.Dimension(neq_displ, neq_plast);
	else
		fNEEMatrixOD.Dimension(neq_plast, neq_displ);
}

/* construct nodal field */
void MFGPElementSupportT::SetNodalField(const dArray2DT& dof)
{
	/* reconstruct displacement field */
	iArrayT nodes;
	fMFShapes->NodalField(dof, fNodalU, nodes);

	/* map range and shift */
	int max;
	nodes.MinMax(fMapShift, max);
	int range = max - fMapShift + 1;
	
	/* dimension */
	fGlobalToNodesUsedMap.Dimension(range);
	fGlobalToNodesUsedMap = -1;

	/* make map */
	for (int i = 0; i < nodes.Length(); i++)
		fGlobalToNodesUsedMap[nodes[i] - fMapShift] = i;
		
	/* set flag */
	fFieldSet = true;
}

void MFGPElementSupportT::GetNodalField(const dArray2DT& dof,
	const iArrayT& nodes, dArray2DT& field) const
{
	/* retrieve stored values */
	if (fFieldSet)
	{
#if __option(extended_errorcheck)
		/* field data might be "free"-ed */
		if (nodes.Length() > fGlobalToNodesUsedMap.Length() ||
		    field.MinorDim() != fNodalU.MinorDim() ||
		    fMapShift < 0) throw ExceptionT::kGeneralFail;
#endif	

		for (int i = 0; i < nodes.Length(); i++)
		{
			int dex = fGlobalToNodesUsedMap[nodes[i] - fMapShift];
			field.SetRow(i, fNodalU(dex));
		}
	}
	/* compute select values */
	else fMFShapes->SelectedNodalField(dof, nodes, field);
}

void MFGPElementSupportT::FreeNodalField(void)
{
	/* unset flag */
	fFieldSet = false;

	/* free memory */
	fNodalU.Free();
	fGlobalToNodesUsedMap.Free();
	fMapShift = -1;
}

/* mark "dead" cells - no active equations */
int MFGPElementSupportT::MarkActiveCells(AutoArrayT<ElementCardT>& elem_cards)
{
	/* loop over cells */
	int active_count = 0;
	for (int i = 0; i < elem_cards.Length(); i++)
	{
		/* equation data */
		int num_eq = fElemEqnosEX.MinorDim(i);
		const int* eq = fElemEqnosEX(i);

		/* look for active equations */
		bool active = false;
		for (int j = 0; j < num_eq && !active; j++)
			if (*eq++ > 0) active = true; //OFFSET

		/* mark cell */
		if (active)
		{
			elem_cards[i].Flag() = 1;
			active_count++;
		}
		else
			elem_cards[i].Flag() = 0;
	}
	return active_count;
}

/* write data for any cell containing the specified node as
* well as the nodes own neighborhood. (map == NULL) means no map. */
void MFGPElementSupportT::TraceNode(ostream& out, int node, const ElementBaseT& element_group)
{
	int d_width = out.precision() + kDoubleExtra;
	out << "\n MFGPElementSupportT::TraceNode: " << node + 1 << endl;

	/* node map */
	const ArrayT<int>* node_map = element_group.ElementSupport().NodeMap();

	/* shape function data */
	MeshFreeSupportT& mf_support = fMFShapes->MeshFreeSupport();

	/* nodal neighborhood information */
	out << " nodal neighborhood:\n";

	/* shape function data */
	iArrayT neighbors;
	dArrayT phi;
	dArray2DT Dphi;
	mf_support.LoadNodalData(node, neighbors, phi, Dphi);

	/* nodal support size */
	dArray2DT nodal_params(neighbors.Length(), mf_support.NodalParameters().MinorDim());
	mf_support.GetSupportParameters(neighbors, nodal_params);

	/* write */
	out << setw(kIntWidth) << "node"
        << setw(d_width) << "phi"
	    << setw(nodal_params.MinorDim()*d_width) << "d_max" << '\n';
	for (int i = 0; i < neighbors.Length(); i++)
	{
		out << setw(kIntWidth) <<
			((node_map != NULL) ? (*node_map)[neighbors[i]] : neighbors[i]) + 1
			<< setw(d_width) << phi[i];
		nodal_params.PrintRow(i, out);
		out << '\n';
	}

	/* integration cell information */
	int num_ip = fMFShapes->NumIP();
	for (int j = 0; j < fUNodeLists.Length(); j++)
	{
		const iArrayT& u_nodes = fUNodeLists[j];
		if (u_nodes.HasValue(node))
		{
			/* cell map */
			const StringT& block_ID = element_group.ElementBlockID(j);
			const iArrayT* element_map = element_group.ElementSupport().ElementMap(block_ID);
		
			out << "    block ID: " << block_ID << '\n';
			out << "(local) cell: " << j + 1 << '\n';

//			out << "     cell: " << ((element_map != NULL) ? (*element_map)[j]: j) + 1 << '\n';
//NOTE: to do the local cell, would need the block size of shift j
	
			/* shape function data */
			iArrayT neighbors;
			dArray2DT phi;
			ArrayT<dArray2DT> Dphi(num_ip);
			mf_support.LoadElementData(j, neighbors, phi, Dphi);
	
			/* nodal support size */
			dArray2DT nodal_params(neighbors.Length(), mf_support.NodalParameters().MinorDim());
			mf_support.GetSupportParameters(neighbors, nodal_params);

			/* write header */
			out << setw(kIntWidth) << "node"
			    << setw(nodal_params.MinorDim()*d_width) << "d_max";
			for (int i = 0; i < num_ip; i++)
				out << setw(d_width - 2) << "phi[" << i+1 << "]";
			out << '\n';
			
			for (int k = 0; k < neighbors.Length(); k++)
			{
				out << setw(kIntWidth) 
				    << ((node_map != NULL) ? (*node_map)[neighbors[k]] : neighbors[k]) + 1;
				nodal_params.PrintRow(k, out);
				for (int i = 0; i < num_ip; i++)
					out << setw(d_width) << phi(i,k);			
				out << '\n';
			}
		}
	}
}

/* weight nodes */
void MFGPElementSupportT::WeightNodes(iArrayT& weight) const
{
	const iArrayT& nodes_used = fMFShapes->MeshFreeSupport().NodesUsed();
	for (int i = 0; i < nodes_used.Length(); i++)
	{
		int& w = weight[nodes_used[i]];
		if (w < 2) w = 2;	
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* class specific data */
void MFGPElementSupportT::CollectNodesData(ostream& out, int max_node_num,
	ModelManagerT* model)
{
	const char caller[] = "MFGPElementSupportT::CollectNodesData";

	/* temp space */
	ArrayT<StringT> set_ID;

	/* meshfree nodes not on the integration grid */
	if (fOffGridID) {
		StringListT::Extract(*fOffGridID, set_ID);
		model->ManyNodeSets(set_ID, fOffGridNodes);
		fOffGridID = NULL;
	} else
		fOffGridNodes.Dimension(0);
	out << "\n Number of nodes off the integration grid. . . . = " << fOffGridNodes.Length() << '\n';

	if (fOffGridNodes.Length() > 0) /* skip empty sets */
	    if (fOffGridNodes.Min() < 0 || fOffGridNodes.Max() >= max_node_num) /* check */
	    	ExceptionT::BadInputValue(caller, "off grid node is out of range: %d > %d",
	    		fOffGridNodes.Max()+1, max_node_num);
	    		
	
	/* interpolant nodes */
	if (fInterpolantID) {
		StringListT::Extract(*fInterpolantID, set_ID);
		model->ManyNodeSets(set_ID, fFENodes);
		fInterpolantID = NULL;
	} else
		fFENodes.Dimension(0);
	out << " Number of interpolant shape function nodes. . . = " << fFENodes.Length() << '\n';

	if (fFENodes.Length() > 0) /* skip empty sets */
	    if (fFENodes.Min() < 0 || fFENodes.Max() >= max_node_num) /* check */
	    	ExceptionT::BadInputValue(caller, "interpolant node is out of range: %d > %d",
	    		fFENodes.Max()+1, max_node_num);


	/* forced meshless nodes */
	if (fMeshlessID) {
		StringListT::Extract(*fMeshlessID, set_ID);
		model->ManyNodeSets(set_ID, fEFGNodes);
		fMeshlessID = NULL;
	} else
		fEFGNodes.Dimension(0);
	out << " Number of pure EFG shape function nodes . . . . = "  << fEFGNodes.Length() << '\n';

	if (fEFGNodes.Length() > 0) /* skip empty sets */
		if (fEFGNodes.Min() < 0 || fEFGNodes.Max() > max_node_num) /* check */
	    	ExceptionT::BadInputValue(caller, "meshless node is out of range: %d > %d",
	    		fEFGNodes.Max()+1, max_node_num);
}

/* set nodes which will have nodally exact shapefunctions */
void MFGPElementSupportT::SetAllFENodes(const iArrayT& fe_nodes)
{
	/* temp collection space */
	iAutoArrayT all_fe_nodes;
	
	/* collect */
	all_fe_nodes.Append(fe_nodes);
	all_fe_nodes.AppendUnique(fFENodes);
	
	/* mark specified meshfree nodes with -1 */
	for (int j = 0; j < fEFGNodes.Length(); j++)
		all_fe_nodes.ChangeValue(fEFGNodes[j], -1);
		//NOTE: this could be more efficient

	/* generate final list */
	fAllFENodes.Dimension(all_fe_nodes.Length() - all_fe_nodes.Count(-1));
	int* from = all_fe_nodes.Pointer();
	int*   to = fAllFENodes.Pointer();
	for (int k = 0; k < all_fe_nodes.Length(); k++)
	{
		if (*from != -1) *to++ = *from;
		from++;
	}
}
