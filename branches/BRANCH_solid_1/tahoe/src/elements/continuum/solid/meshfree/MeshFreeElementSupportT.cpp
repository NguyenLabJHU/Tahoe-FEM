/* $Id: MeshFreeElementSupportT.cpp,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (11/12/1999)                                          */

#include "MeshFreeElementSupportT.h"

#include "fstreamT.h"
#include "iAutoArrayT.h"
#include "MeshFreeShapeFunctionT.h"
#include "ElementCardT.h"
#include "MeshFreeSupportT.h"

#include "ModelFileT.h"
#include "ExodusT.h"

/* needed for TraceNode */
#include "FEManagerT.h"
#include "ElementBaseT.h"

/* parameters */
const int kHeadRoom = 10; // percent

#ifdef __MPI__
#include "mpi.h"
#endif

/* constructor */
MeshFreeElementSupportT::MeshFreeElementSupportT(ifstreamT& in):
	fMFShapes(NULL),
	fLocGroup(kHeadRoom),
	fNumElemenNodes(0),
	fNEEArray(kHeadRoom),
	fNEEMatrix(kHeadRoom),
	fFieldSet(false),
	fMapShift(-1)
{
	/* read */
	in >> fMeshFreeCode;
	in >> fd_max;
	in >> fComplete;
	in >> fStoreShape;
	in >> fAutoBorder;

	/* check values */
	if (fd_max < 1.0) throw eBadInputValue;
	if (fComplete < 1) throw eBadInputValue;
	if (fStoreShape != 0 && fStoreShape != 1) throw eBadInputValue;
	if (fAutoBorder != 0 && fAutoBorder != 1) throw eBadInputValue;

#ifdef __MPI__
	//TEMP
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size > 1 && fAutoBorder)
		cout << "\n ::MeshFreeElementSupportT: AutoBorder not extended to parallel" << endl;
#endif
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
	out << " Mesh-free formulation . . . . . . . . . . . . . = " << fMeshFreeCode << '\n';
	out << "    eq. " << MeshFreeSupportT::kEFG  << ", EFG\n";
	out << "    eq. " << MeshFreeSupportT::kRKPM << ", RKPM\n";
	out << " Domain of influence scale factor (d_max). . . . = " << fd_max        << '\n';
	out << " Order of completeness of basis functions. . . . = " << fComplete     << '\n';
	out << " Store all shape functions and derivatives . . . = " << fStoreShape   << '\n';
	out << " Auto selection/generation of border transition. = " << fAutoBorder;
	out << ((fAutoBorder == 1) ? " (ACTIVE)" : " (INACTIVE)") << '\n';
	out.flush();
}

/* initialization */
void MeshFreeElementSupportT::InitSupport(ifstreamT& in, ostream& out,
	AutoArrayT<ElementCardT>& elem_cards, const iArrayT& surface_nodes,
	int numDOF, int max_node_num, const StringT& model_file,
	IOBaseT::FileTypeT format)
{
	/* configure variable length element arrays */
	fElemNodesEX = &(fMFShapes->ElementNeighbors());
	fElemEqnosEX.Configure(fMFShapes->ElementNeighborsCounts(), numDOF);
	
	/* set element card pointers */
	int num_cells = elem_cards.Length();
	fUNodeLists.Allocate(num_cells);
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
	EchoNodesData(in, out, max_node_num, model_file, format);

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
	fGlobalToNodesUsedMap.Allocate(range);
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
		    fMapShift < 0) throw eGeneralFail;
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
	const iArrayT* node_map = element_group.FEManager().NodeMap();

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
	dArrayT d_max(neighbors.Length());
	mf_support.GetDmax(neighbors, d_max);

	/* write */
	out << setw(kIntWidth) << "node"
	    << setw(d_width) << "d_max"
	    << setw(d_width) << "phi"   << '\n';
	for (int i = 0; i < neighbors.Length(); i++)
	{
		out << setw(kIntWidth) <<
			((node_map != NULL) ? (*node_map)[neighbors[i]] : neighbors[i]) + 1
		    << setw(d_width) << d_max[i]
		    << setw(d_width) << phi[i] << '\n';
	}

	/* integration cell information */
	int num_ip = fMFShapes->NumIP();
	for (int j = 0; j < fUNodeLists.Length(); j++)
	{
		const iArrayT& u_nodes = fUNodeLists[j];
		if (u_nodes.HasValue(node))
		{
			/* cell map */
			int block_ID = element_group.ElementBlockID(j);
			const iArrayT* element_map = element_group.FEManager().ElementMap(block_ID);
		
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
			dArrayT d_max(neighbors.Length());
			mf_support.GetDmax(neighbors, d_max);

			/* write header */
			out << setw(kIntWidth) << "node"
			    << setw(d_width) << "d_max";
			for (int i = 0; i < num_ip; i++)
				out << setw(d_width - 2) << "phi[" << i+1 << "]";
			out << '\n';
			
			for (int k = 0; k < neighbors.Length(); k++)
			{
				out << setw(kIntWidth) <<
					((node_map != NULL) ? (*node_map)[neighbors[k]] : neighbors[k]) + 1
					<< setw(d_width) << d_max[k];
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
	const StringT& model_file, IOBaseT::FileTypeT format)
{
	/* could echo "sampling points" */

	/* EFG nodes not on the integration grid */
	int num_offgrid;
	in >> num_offgrid;
	out << "\n Number of nodes off the integration grid. . . . = "
	    <<num_offgrid << '\n';
	if (num_offgrid > 0)
	{
		/* read data */
		ReadNodesData(in, num_offgrid, model_file, format, fOffGridNodes);

		/* correct offset */
		if (fOffGridNodes.Length() > 0) // empty sets
		{
			fOffGridNodes--;

			/* check */
			if (fOffGridNodes.Min() < 0 || fOffGridNodes.Max() >= max_node_num)
			{
				cout << "\n MeshFreeElementSupportT::EchoNodesData: off grid EFG node out of range" << endl;
				throw eBadInputValue;
			}
		}
	}
	
	/* interpolant nodes */
	int num_FE;
	in >> num_FE;
	out << " Number of interpolant shape function nodes. . . = " << num_FE << '\n';
	if (num_FE > 0)
	{
		/* read data */
		ReadNodesData(in, num_FE, model_file, format, fFENodes);

		/* correct offset */
		if (fFENodes.Length() > 0) // empty sets
		{
			fFENodes--;

			/* check */
			if (fFENodes.Min() < 0 || fFENodes.Max() >= max_node_num)
			{
				cout << "\n MeshFreeElementSupportT::EchoNodesData: interpolant node out of range" << endl;
				throw eBadInputValue;
			}
		}
	}

	/* pure EFG nodes */
	int num_EFG;
	in >> num_EFG;
	out << " Number of pure EFG shape function nodes . . . . = " << num_EFG << '\n';
	if (num_EFG > 0)
	{
		/* read data */
		ReadNodesData(in, num_EFG, model_file, format, fEFGNodes);

		/* correct offset */
		if (fEFGNodes.Length() > 0) // empty sets
		{
			fEFGNodes--;

			/* check */
			if (fEFGNodes.Min() < 0 || fEFGNodes.Max() > max_node_num)
			{
				cout << "\n MeshFreeElementSupportT::EchoNodesData: set EFG node out of range" << endl;
				throw eBadInputValue;
			}
		}
	}
}

void MeshFreeElementSupportT::ReadNodesData(ifstreamT& in, int num_id, const StringT& model_file,
	IOBaseT::FileTypeT format, iArrayT& nodes)
{
	/* read contact nodes */
	switch (format)
	{
		case IOBaseT::kTahoe:
		{
			nodes.Allocate(num_id);
			in >> nodes;
			break;
		}
		case IOBaseT::kTahoeII:
		{
			/* number of node sets */
			if (num_id > 0)
			{
				/* open database */
				ModelFileT database;
				database.OpenRead(model_file);

				/* echo set ID's */
				iArrayT ID_list(num_id);
				in >> ID_list;
				
				/* collect */
				if (database.GetNodeSets(ID_list, nodes) != ModelFileT::kOK)
					throw eBadInputValue;
			}
			break;
		}
		case IOBaseT::kExodusII:
		{
			/* number of node sets */
			if (num_id > 0)
			{
				/* echo set ID's */
				iArrayT ID_list(num_id);
				in >> ID_list;

				/* open database */
				ExodusT database(cout);
				database.OpenRead(model_file);
				
				/* read collect all nodes in sets */
				database.ReadNodeSets(ID_list, nodes);
			}				
			break;
		}
		default:

			cout << "\n MeshFreeElementSupportT::ReadNodesData: unsupported input format: ";
			cout << format << endl;
			throw eGeneralFail;
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
	fAllFENodes.Allocate(all_fe_nodes.Length() - all_fe_nodes.Count(-1));
	int* from = all_fe_nodes.Pointer();
	int*   to = fAllFENodes.Pointer();
	for (int k = 0; k < all_fe_nodes.Length(); k++)
	{
		if (*from != -1) *to++ = *from;
		from++;
	}
}
