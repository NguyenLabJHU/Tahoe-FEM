/* $Id: MeshFreeElementSupportT.cpp,v 1.14 2004-06-17 07:41:26 paklein Exp $ */
/* created: paklein (11/12/1999) */
#include "MeshFreeElementSupportT.h"

#include "ifstreamT.h"
#include "iAutoArrayT.h"
#include "MeshFreeShapeFunctionT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "ElementCardT.h"
#include "MeshFreeSupportT.h"
#include "ElementBaseT.h"

#include "ModelManagerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* parameters */
const int kHeadRoom = 10; // percent

/* constructor */
MeshFreeElementSupportT::MeshFreeElementSupportT(ifstreamT& in):
	fMFShapes(NULL),
	fNodalShapes(NULL),
	fLocGroup(kHeadRoom),
	fNumElemenNodes(0),
	fNEEArray(kHeadRoom, true),
	fNEEMatrix(kHeadRoom, true),
	fFieldSet(false),
	fMapShift(-1)
{
	/* read */
	in >> fAutoBorder;

	/* check values */
	if (fAutoBorder != 0 && fAutoBorder != 1) throw ExceptionT::kBadInputValue;
}

/* accessors */
MeshFreeSupportT& MeshFreeElementSupportT::MeshFreeSupport(void) const
{
	return fMFShapes->MeshFreeSupport();
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void MeshFreeElementSupportT::PrintControlData(ostream& out) const
{
	/* echo */
	out << " Auto selection/generation of border transition. = " << fAutoBorder;
	out << ((fAutoBorder == 1) ? " (ACTIVE)" : " (INACTIVE)") << '\n';
	out.flush();
}

/* initialization */
void MeshFreeElementSupportT::InitSupport(ifstreamT& in, ostream& out,
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
	
	/* echo FE/meshfree nodes */
	EchoNodesData(in, out, max_node_num, model);

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
int MeshFreeElementSupportT::SetElementNodes(int element)
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

/* construct nodal field */
void MeshFreeElementSupportT::SetNodalField(const dArray2DT& dof)
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

void MeshFreeElementSupportT::GetNodalField(const dArray2DT& dof,
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

void MeshFreeElementSupportT::FreeNodalField(void)
{
	/* unset flag */
	fFieldSet = false;

	/* free memory */
	fNodalU.Free();
	fGlobalToNodesUsedMap.Free();
	fMapShift = -1;
}

/* mark "dead" cells - no active equations */
int MeshFreeElementSupportT::MarkActiveCells(AutoArrayT<ElementCardT>& elem_cards)
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
void MeshFreeElementSupportT::TraceNode(ostream& out, int node, const ElementBaseT& element_group)
{
	int d_width = out.precision() + kDoubleExtra;
	out << "\n MeshFreeElementSupportT::TraceNode: " << node + 1 << endl;

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
void MeshFreeElementSupportT::WeightNodes(iArrayT& weight) const
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
void MeshFreeElementSupportT::EchoNodesData(ifstreamT& in, ostream& out, int max_node_num,
	ModelManagerT* model)
{
	/* could echo "sampling points" */

	/* temp space */
	ArrayT<StringT> set_ID;

	/* EFG nodes not on the integration grid */
	model->NodeSetList (in, set_ID);
	model->ManyNodeSets (set_ID, fOffGridNodes);
	out << "\n Number of nodes off the integration grid. . . . = "
	    << fOffGridNodes.Length() << '\n';

	if (fOffGridNodes.Length() > 0) // empty sets
	  {
	    /* check */
	    if (fOffGridNodes.Min() < 0 || fOffGridNodes.Max() >= max_node_num)
	      {
		cout << "\n MeshFreeElementSupportT::EchoNodesData: off grid EFG node out of range" << endl;
		throw ExceptionT::kBadInputValue;
	      }
	  }
	
	/* interpolant nodes */
	model->NodeSetList (in, set_ID);
	model->ManyNodeSets (set_ID, fFENodes);
	out << " Number of interpolant shape function nodes. . . = " 
	    << fFENodes.Length() << '\n';

	if (fFENodes.Length() > 0) // empty sets
	  {
	    /* check */
	    if (fFENodes.Min() < 0 || fFENodes.Max() >= max_node_num)
	      {
		cout << "\n MeshFreeElementSupportT::EchoNodesData: interpolant node out of range" << endl;
		throw ExceptionT::kBadInputValue;
	      }
	  }

	/* pure EFG nodes */
	model->NodeSetList (in, set_ID);
	model->ManyNodeSets (set_ID, fEFGNodes);
	out << " Number of pure EFG shape function nodes . . . . = " 
	    << fEFGNodes.Length() << '\n';

	if (fEFGNodes.Length() > 0) // empty sets
	  {
	    /* check */
	    if (fEFGNodes.Min() < 0 || fEFGNodes.Max() > max_node_num)
	      {
		cout << "\n MeshFreeElementSupportT::EchoNodesData: set EFG node out of range" << endl;
		throw ExceptionT::kBadInputValue;
	      }
	  }
}

/* set nodes which will have nodally exact shapefunctions */
void MeshFreeElementSupportT::SetAllFENodes(const iArrayT& fe_nodes)
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
