/* $Id: ElementBaseT.cpp,v 1.7.2.2 2001-10-15 19:38:30 sawimme Exp $ */
/* created: paklein (05/24/1996)                                          */

#include "ElementBaseT.h"

#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"

/* array behavior */
const bool ArrayT<const RaggedArray2DT<int>*>::fByteCopy = true;

/* constructor */
ElementBaseT::ElementBaseT(FEManagerT& fe_manager):
	fFEManager(fe_manager),
	fNodes(NULL),
	fController(NULL),
	fNumElemNodes(0),
	fNumElemEqnos(0),
	fNumElements(0),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{
	/* get pointer to the node manager */
	fNodes  = fFEManager.NodeManager();
	fNumSD  = fNodes->NumSD();
	fNumDOF = fNodes->NumDOF(); // NOTE: won't be good for multi-physics

	fAnalysisCode = fFEManager.Analysis();
}

/* destructor */
ElementBaseT::~ElementBaseT(void) {	}

/* run status */
const GlobalT::StateT& ElementBaseT::RunState(void) const
{ return fFEManager.RunState(); }

/* allocates space and reads connectivity data */
void ElementBaseT::Initialize(void)
{
	/* set console variables */
	int index = fFEManager.ElementGroupNumber(this) + 1;
	StringT name;
	name.Append(index);
	name.Append("_element_group");
	iSetName(name);
	iAddVariable("num_elements", *((const int*) &fNumElements));
	iAddVariable("num_element_nodes", *((const int*) &fNumElemNodes));

	/* streams */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();

	/* control data */
	PrintControlData(out);

	/* element connectivity data */
	EchoConnectivityData(in, out);

	/* dimension */
	fLHS.Allocate(fNumElemEqnos);	
	fRHS.Allocate(fNumElemEqnos);
}

/*
* Re-initialize: signal to element group that the global
* equations numbers are going to be reset so that the group
* has the opportunity to reconnect and should reinitialize
* an dependencies on global equation numbers obtained from the
* NodeManagerT.
*
* NOTE: any memory allocated after initial construction (or since
* the last Reinitialize) should be "shuffled down" at this point, ie.
* reallocated and copied, to make room for the global stiffness
* matrix.
*/
void ElementBaseT::Reinitialize(void)
{
	/* do nothing by default */
}

/* set the controller */
void ElementBaseT::SetController(eControllerT* controller)
{
	fController = controller;
}

/* form of tangent matrix - symmetric by default */
GlobalT::SystemTypeT ElementBaseT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* solution calls */
void ElementBaseT::FormLHS(void)
{
	try { LHSDriver(); }
	catch (int error)
	{
		cout << "\n ElementBaseT::FormLHS: " << fFEManager.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fFEManager.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{
			ostream& out = fFEManager.Output();
		
			/* header */
			out << "\n ElementBaseT::FormLHS: caught exception " << error << '\n';
			out <<   "      Time: " << fFEManager.Time() << '\n';
			out <<   "      Step: " << fFEManager.StepNumber() << '\n';
			out <<   " Time step: " << fFEManager.TimeStep() << '\n';
		
			/* write current element information to main out */
			CurrElementInfo(out);
			cout << "     See output file for current element information\n";
		}
		else
			cout << "     Current element information not available\n";
		cout.flush();
		throw error;
	}
}

void ElementBaseT::FormRHS(void)
{
	try { RHSDriver(); }
	catch (int error)
	{
		cout << "\n ElementBaseT::FormRHS: " << fFEManager.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fFEManager.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{
			ostream& out = fFEManager.Output();
		
			/* header */
			out << "\n ElementBaseT::FormRHS: caught exception " << error << '\n';
			out <<   "      Time: " << fFEManager.Time() << '\n';
			out <<   "      Step: " << fFEManager.StepNumber() << '\n';
			out <<   " Time step: " << fFEManager.TimeStep() << '\n';
		
			/* write current element information to main out */
			CurrElementInfo(out);
			cout << "     See output file for current element information\n";
		}
		else
			cout << "     Current element information not available\n";
		cout.flush();
		throw error;
	}
}

/* initialize/finalize time increment */
void ElementBaseT::InitStep(void) { }
void ElementBaseT::CloseStep(void) { }

/* resets to the last converged solution */
void ElementBaseT::ResetStep(void)
{
	/* do nothing by default */
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT ElementBaseT::RelaxSystem(void)
{
	/* no relaxation */
	return GlobalT::kNoRelax;
}

/* append element equations numbers to the list */
void ElementBaseT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* get local equations numbers */
	fNodes->SetLocalEqnos(fConnectivities, fEqnos);

	/* add to list */
	eq_1.Append(&fEqnos);
}

/* appends group connectivities to the array (X -> geometry, U -> field) */
void ElementBaseT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
	/* append connectivities */
	connects.AppendUnique(&fConnectivities);
}

void ElementBaseT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_2)
	/* by default field nodes are geometry nodes */
	ConnectsX(connects_1);
}

/* returns a pointer to the specified LoadTime function */
LoadTime* ElementBaseT::GetLTfPtr(int num) const
{
	return fFEManager.GetLTfPtr(num);
}

/* initial condition/restart functions
*
* Set to initial conditions.  The restart functions
* should read/write any data that overrides the default
* values */
void ElementBaseT::InitialCondition(void)
{
	//do nothing
}

void ElementBaseT::ReadRestart(istream& in)
{
	/* stream check */
	if (!in.good()) throw eGeneralFail;
}

void ElementBaseT::WriteRestart(ostream& out) const
{
	/* stream check */
	if (!out.good()) throw eGeneralFail;
}

/* returns 1 if DOF's are interpolants of the nodal values */
int ElementBaseT::InterpolantDOFs(void) const { return 1; }

void ElementBaseT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
#if __option(extended_errorcheck)
	if (nodes.Length() != DOFs.MajorDim() ||
	    DOFs.MinorDim() != fNumDOF) throw eSizeMismatch;

#endif

	const dArray2DT& all_DOFs = fNodes->Displacements();
	DOFs.RowCollect(nodes, all_DOFs); // no check of nodes used

//NOTE - This function is added only for completeness. If the
//       DOF's are interpolant, there should be no reason to
//       collect the displacements in this way. Use SetLocalU
}

/* block ID for the specified element */
int ElementBaseT::ElementBlockID(int element) const
{
	if (element < 0 || element >= fNumElements) throw eOutOfRange;
	
	int blockID = 0;
	bool found = false;
	for (int i = 0; i < fBlockData.MajorDim() && !found; i++)
		if (element >= fBlockData(i, kStartNum) &&
		    element <  fBlockData(i, kStartNum) + fBlockData(i, kBlockDim))
		{
			blockID = fBlockData(i, kID);
			found = true;
		}
	if (!found) throw eGeneralFail;
	return blockID;
}

/* weight the computational effort of every node */
void ElementBaseT::WeightNodalCost(iArrayT& weight) const
{
	int base_weight = 1;
	int* p = fConnectivities.Pointer();
	int n = fConnectivities.Length();
	for (int i = 0; i < n; i++)
	{
		if (weight[*p] < base_weight) weight[*p] = base_weight;
		p++;
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* get local element data, X for geometry, U for
* field variables */
const LocalArrayT& ElementBaseT::SetLocalX(LocalArrayT& localarray)
{
	localarray.SetLocal(CurrentElement().NodesX());
	return localarray;
}

const LocalArrayT& ElementBaseT::SetLocalU(LocalArrayT& localarray)
{
	localarray.SetLocal(CurrentElement().NodesU());
	return localarray;
}

/* assembling the left and right hand sides */
void ElementBaseT::AssembleRHS(void) const
{
	fFEManager.AssembleRHS(fRHS, CurrentElement().Equations());
}

void ElementBaseT::AssembleLHS(void) const
{
	fFEManager.AssembleLHS(fLHS, CurrentElement().Equations());
}

/* print element group data */
void ElementBaseT::PrintControlData(ostream& out) const
{
#pragma unused(out)
}

/* echo element connectivity data, resolve material pointers
* and set the local equation numbers */
void ElementBaseT::EchoConnectivityData(ifstreamT& in, ostream& out)
{	
	out << "\n Element Connectivity:\n";
	
	/* read */
	ReadConnectivity(in, out);

	/* write */
	WriteConnectivity(out);
}

/* resolve input format types */
void ElementBaseT::ReadConnectivity(ifstreamT& in, ostream& out)
{	
	/* read from parameter file */
        iArrayT indexes, matnums;
	ModelManagerT* model = fFEManager.ModelManager();
	model->ElementBlockList (in, indexes, matnums);

	/* allocate block map */
	int num_blocks = indexes.Length();
	fBlockData.Allocate(num_blocks, kBlockDataSize);

	/* read from parameter file */
	int elem_count = 0;
	int nen;
	for (int b=0; b < num_blocks; b++)
	  {
	    /* check number of nodes */
	    int num_elems, num_nodes;
	    model->ElementGroupDimensions (indexes[b], num_elems, num_nodes);
	    if (b == 0)
	      nen = num_nodes;
	    else if (nen != num_nodes)
	      {
		cout << "\n ElementBaseT::ReadConnectivity: minor dimension "
		     << num_nodes << " of block " << b+1 << '\n';
		cout <<   "     does not match dimension of previous blocks "
		     << nen << endl;
		throw eBadInputValue;
	      }
	    
	    /* store block data */
	    fBlockData (b, kID) = indexes[b] + 1; // use global index as ID value
	    fBlockData (b, kStartNum) = elem_count;
	    fBlockData (b, kBlockDim) = num_elems;
	    fBlockData (b, kBlockMat) = matnums[b] - 1; // offset

	    /* increment element count */
	    elem_count += num_elems;
	  }
	
	/* write to single list */
	if (num_blocks == 1)
	  fConnectivities = model->ElementGroup (fBlockData (0, kID) - 1);
	else
	  {
	    /* allocate */
	    fConnectivities.Allocate (elem_count, nen);

	    /* by-block */
	    for (int i=0; i < fBlockData.MajorDim(); i++)
	      fConnectivities.BlockRowCopyAt (model->ElementGroup (indexes[i]), fBlockData (i, kStartNum));
	  }

	/* set dimensions */
	fNumElements  = fConnectivities.MajorDim();
	fNumElemNodes = fConnectivities.MinorDim();
	
	/* connectivity returned empty */
	if (fNumElemNodes == 0)
	{
		fNumElemNodes = DefaultNumElemNodes();
		fConnectivities.Allocate(0, fNumElemNodes);
	}

	/* derived dimensions */	
	fNumElemEqnos = fNumElemNodes*fNumDOF;
	fEqnos.Allocate(fNumElements, fNumElemEqnos);

	/* set pointers in element cards */
	SetElementCards();
}

/* resolve output formats */
void ElementBaseT::WriteConnectivity(ostream& out) const
{	
	out << " Number of elements. . . . . . . . . . . . . . . = " << fNumElements  << '\n';

	/* verbose output */
	if (fFEManager.PrintInput())
	{
		/* write header */
		out << setw(kIntWidth) << "no.";
		out << setw(kIntWidth) << "mat.";
		for (int j = 1; j <= ((fNumElemNodes < 9) ? fNumElemNodes : 8); j++)
		{
			int numwidth = (j < 10) ? 1 : ((j < 100) ? 2 : 3);		
			out << setw(kIntWidth - (numwidth + 1)) << "n[";
			out << j << "]";
		}
		out << endl;
				
		/* write material number and connectivity */
		iArrayT nodesX(fNumElemNodes);
		for (int i = 0; i < fNumElements; i++)
		{
			const ElementCardT& elcard = fElementCards[i];
		
			out << setw(kIntWidth) << i + 1;
			out << setw(kIntWidth) << elcard.MaterialNumber() + 1;		
		
			/* nodes defining the geometry */
			nodesX = elcard.NodesX();
			nodesX++;
			out << nodesX.wrap(8, kIntWidth) << '\n';
			nodesX--;
		}
		out << endl;
	}
}

/* generate connectivities with local numbering -
* returns the number of nodes used by the element group */
int ElementBaseT::MakeLocalConnects(iArray2DT& localconnects)
{
	/* compressed number range */
	int min   = fConnectivities.Min();
	int range = fConnectivities.Max() - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int i = 0; i < fConnectivities.Length(); i++)
		node_map[fConnectivities[i] - min] = 1;

	/* set node map */
	int localnum = 0;
	for (int j = 0; j < node_map.Length(); j++)
	    if (node_map[j] == 1)
		    node_map[j] = localnum++;

	/* connectivities with local node numbering */
localconnects.Allocate(fNumElements, fNumElemNodes);
for (int k = 0; k < localconnects.Length(); k++)
	   localconnects[k] = node_map[fConnectivities[k] - min];

	return localnum;
}

void ElementBaseT::NodesUsed(ArrayT<int>& nodes_used) const
{
	/* compressed number range */
	int min   = fConnectivities.Min();
	int range = fConnectivities.Max() - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int i = 0; i < fConnectivities.Length(); i++)
		node_map[fConnectivities[i] - min] = 1;

	/* collect list */
	nodes_used.Allocate(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodes_used[dex++] = j + min;
}

/* return pointer to block data given the ID */
const int* ElementBaseT::BlockData(int block_ID) const
{
	/* resolve block ID */
	int block_num = -1;
	for (int j = 0; j < fBlockData.MajorDim() && block_num == -1; j++)
		if (fBlockData(j, kID) == block_ID) block_num = j;

	/* check */
	if (block_num == -1)
	{
		cout << "\n ElementBaseT::BlockData: block ID ";
		cout << block_ID << " not found in\n";
		cout <<   "     element group " << fFEManager.ElementGroupNumber(this) + 1;
		cout << ". Block data:\n";
		cout << setw(kIntWidth) << "ID"
		     << setw(kIntWidth) << "start"
		     << setw(kIntWidth) << "size"
		     << setw(kIntWidth) << "mat." << '\n';
		cout << fBlockData << endl;
		throw eBadInputValue;
	}

	/* return */
	return fBlockData(block_num);
}

/* write all current element information to the stream */
void ElementBaseT::CurrElementInfo(ostream& out) const
{
	if (!fElementCards.InRange()) return;
	
	out << "\n   element group: " << fFEManager.ElementGroupNumber(this) + 1 << '\n';
	out << "\n current element: " << fElementCards.Position() + 1 << '\n';
	
	/* node number map */
	const iArrayT* node_map = fFEManager.NodeMap();
	iArrayT temp;

	out <<   " connectivity(x):\n";
	if (node_map)
	{
		const iArrayT& nodes_X = CurrentElement().NodesX();
		temp.Allocate(nodes_X.Length());
		for (int i = 0; i < nodes_X.Length(); i++)
			temp[i] = (*node_map)[nodes_X[i]];
	}
	else
		temp = CurrentElement().NodesX();
	temp++;
	out << temp.wrap(4) << '\n';

	out <<   " connectivity(u):\n";
	if (node_map)
	{
		const iArrayT& nodes_U = CurrentElement().NodesU();
		temp.Allocate(nodes_U.Length());
		for (int i = 0; i < nodes_U.Length(); i++)
			temp[i] = (*node_map)[nodes_U[i]];
	}
	else
		temp = CurrentElement().NodesU();
	temp++;
	out << temp.wrap(4) << '\n';

	out <<   " equations:\n";
	out << (CurrentElement().Equations()).wrap(4) << '\n';
}

/* set element cards array */
void ElementBaseT::SetElementCards(void)
{
	/* allocate */
	fElementCards.Allocate(fNumElements);

	/* loop over blocks to set pointers */
	int numberofnodes = fNodes->NumNodes();
	int count = 0;
	for (int i = 0; i < fBlockData.MajorDim(); i++)
	{
		int dim = fBlockData(i, kBlockDim);
		int mat = fBlockData(i, kBlockMat);

		for (int j = 0; j < dim; j++)
		{
			ElementCardT& element_card = fElementCards[count];
	
			/* material number */
			element_card.SetMaterialNumber(mat);

			/* set pointers */
			iArrayT& nodes = element_card.NodesX();
			fConnectivities.RowAlias(count, nodes);
			fEqnos.RowAlias(count, element_card.Equations());
			
			/* check node numbers */
			int min, max;
			nodes.MinMax(min, max);
			if (min < 0 || max >= numberofnodes)
			{
				cout << "\n ElementBaseT::SetElementCards: nodes {" << min + 1
				     << "," << max + 1 << "} in element " << count + 1 << "\n";
				cout <<   "     (" << j + 1 << " in block " <<  i + 1 << ") of group "
				     << fFEManager.ElementGroupNumber(this) + 1 << " are out of range" << endl;
				throw eBadInputValue;
			}		
	
			/* next */
			count++;
		}
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* return the default number of element nodes */
int ElementBaseT::DefaultNumElemNodes(void) const { return 0; }
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.
