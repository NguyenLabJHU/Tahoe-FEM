/* $Id: ElementBaseT.cpp,v 1.12 2002-02-11 01:39:02 paklein Exp $ */
/* created: paklein (05/24/1996) */

#include "ElementBaseT.h"

#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "ModelManagerT.h"
#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "LocalArrayT.h"

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

#if __option(extended_errorcheck)
	if (fConnectivities.Length() != fEqnos.Length()) throw eSizeMismatch;
#endif

	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* get local equations numbers */
		fNodes->SetLocalEqnos(*fConnectivities[i], fEqnos[i]);

		/* add to list of equation numbers */
		eq_1.Append(&fEqnos[i]);
	}
}

/* appends group connectivities to the array (X -> geometry, U -> field) */
void ElementBaseT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
	/* loop over connectivity blocks */
	for (int i = 0; i < fConnectivities.Length(); i++)
		connects.AppendUnique(fConnectivities[i]);
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
const StringT& ElementBaseT::ElementBlockID(int element) const
{
	if (element < 0 || element >= fNumElements) {
		cout << "\n ElementBaseT::ElementBlockID: element number " << element << " is out of range {0,"
		    << fNumElements - 1 << "}" << endl;
		throw eOutOfRange;
	}
	
	bool found = false;
	for (int i = 0; i < fBlockData.Length(); i++)
		if (element >= fBlockData[i].StartNumber() &&
		    element <  fBlockData[i].StartNumber() + fBlockData[i].Dimension())
			return fBlockData[i].ID();

	if (!found) {
		cout << "\n ElementBaseT::ElementBlockID: could not resolve block ID for element "
		     << element << endl;
		throw eGeneralFail;
	}
	return fBlockData[0].ID(); /* dummy */
}

/* weight the computational effort of every node */
void ElementBaseT::WeightNodalCost(iArrayT& weight) const
{
	int base_weight = 1;

	for (int i=0; i < fNumElements; i++)
	  {
	    const iArrayT& elemnodes = fElementCards[i].NodesX();
	    int* p = elemnodes.Pointer();
	    for (int n=0; n < elemnodes.Length(); n++)
	      if (weight[*p] < base_weight) 
		weight[*p] = base_weight;
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
#pragma unused(out)

	/* read from parameter file */
	ArrayT<StringT> elem_ID;
	iArrayT matnums;
	ModelManagerT* model = fFEManager.ModelManager();
	model->ElementBlockList(in, elem_ID, matnums);

	/* allocate block map */
	int num_blocks = elem_ID.Length();
	fBlockData.Allocate(num_blocks);
	fConnectivities.Allocate (num_blocks);

	/* read from parameter file */
	int elem_count = 0;
	int nen = 0;
	for (int b=0; b < num_blocks; b++)
	{
	    /* check number of nodes */
	    int num_elems, num_nodes;
	    model->ElementGroupDimensions(elem_ID[b], num_elems, num_nodes);
	    
	    /* set if unset */
	    if (nen == 0) nen = num_nodes;
	    
	    /* consistency check */
	    if (num_nodes != 0 && nen != num_nodes)
		{
			cout << "\n ElementBaseT::ReadConnectivity: minor dimension "
                 << num_nodes << " of block " << b+1 << '\n';
			cout <<   "     does not match dimension of previous blocks "
                 << nen << endl;
			throw eBadInputValue;
		}
	    
	    /* store block data */
	    fBlockData[b].Set(elem_ID[b], elem_count, num_elems, matnums[b] - 1); // offset

	    /* increment element count */
	    elem_count += num_elems;

	    /* load connectivity from database into model manager */
	    model->ReadConnectivity(elem_ID[b]);

	    /* set pointer to connectivity list */
	    fConnectivities[b] = model->ElementGroupPointer(elem_ID[b]);
	}
	  
	/* set dimensions */
	fNumElements  = elem_count;
	fNumElemNodes = nen;
	
	/* connectivity returned empty */
	if (fNumElemNodes == 0) fNumElemNodes = DefaultNumElemNodes();

	/* derived dimensions */	
	fNumElemEqnos = fNumElemNodes*fNumDOF;
	fEqnos.Allocate (num_blocks);
	for (int be=0; be < num_blocks; be++)
	  {
	    int numblockelems = fConnectivities[be]->MajorDim();
	    fEqnos[be].Allocate(numblockelems, fNumElemEqnos);
	  }

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
       int num_blocks = fBlockData.Length();

       iArrayT mins (num_blocks);
       iArrayT maxes (num_blocks);
       for (int i=0; i < fBlockData.Length(); i++)
	 {
	   mins[i] = fConnectivities[i]->Min();
	   maxes[i] = fConnectivities[i]->Max();
	 }

	/* compressed number range */
	int min   = mins.Min();
	int range = maxes.Max() - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int b=0; b < num_blocks; b++)
	  {
	    const iArray2DT* conn = fConnectivities[b];
	    int *pc = conn->Pointer();
	    for (int i = 0; i < conn->Length(); i++)
	      node_map[*pc++ - min] = 1;
	  }

	/* set node map */
	int localnum = 0;
	for (int j = 0; j < node_map.Length(); j++)
	    if (node_map[j] == 1)
		    node_map[j] = localnum++;

	/* connectivities with local node numbering */
	localconnects.Allocate(fNumElements, fNumElemNodes);
	int *plocal = localconnects.Pointer();
	for (int b=0; b < num_blocks; b++)
	  {
	    const iArray2DT* conn = fConnectivities[b];
	    int *pc = conn->Pointer();
	    for (int i = 0; i < conn->Length(); i++)
	      *plocal++ = node_map [*pc++ - min];
	  }

	return localnum;
}

void ElementBaseT::NodesUsed(ArrayT<int>& nodes_used) const
{
       int num_blocks = fBlockData.Length();

       iArrayT mins (num_blocks);
       iArrayT maxes (num_blocks);
       for (int i = 0; i < fBlockData.Length(); i++)
	 {
	   mins[i] = fConnectivities[i]->Min();
	   maxes[i] = fConnectivities[i]->Max();
	 }

	/* compressed number range */
	int min   = mins.Min();
	int range = maxes.Max() - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int b=0; b < num_blocks; b++)
	  {
	    const iArray2DT* conn = fConnectivities[b];
	    int *pc = conn->Pointer();
	    for (int i = 0; i < conn->Length(); i++)
	      node_map[*pc++ - min] = 1;
	  }

	/* collect list */
	nodes_used.Allocate(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodes_used[dex++] = j + min;
}

/* return pointer to block data given the ID */
const ElementBlockDataT& ElementBaseT::BlockData(const StringT& block_ID) const
{
	/* resolve block ID */
	int block_num = -1;
	for (int j = 0; j < fBlockData.Length() && block_num == -1; j++)
		if (fBlockData[j].ID() == block_ID) block_num = j;

	/* check */
	if (block_num == -1)
	{
		cout << "\n ElementBaseT::BlockData: block ID ";
		cout << block_ID << " not found in\n";
		cout <<   "     element group " << fFEManager.ElementGroupNumber(this) + 1;
		cout << ". Block data:\n";
		cout << setw(12) << "ID"
		     << setw(kIntWidth) << "start"
		     << setw(kIntWidth) << "size"
		     << setw(kIntWidth) << "mat." << '\n';

		for (int i = 0; i < fBlockData.Length(); i++)
			cout << setw(12) << fBlockData[i].ID()
                 << setw(kIntWidth) << fBlockData[i].StartNumber()
                 << setw(kIntWidth) << fBlockData[i].Dimension()
                 << setw(kIntWidth) << fBlockData[i].MaterialID() << '\n';

		cout.flush();
		throw eBadInputValue;
	}

	/* return */
	return fBlockData[block_num];
}

/* write all current element information to the stream */
void ElementBaseT::CurrElementInfo(ostream& out) const
{
	if (!fElementCards.InRange()) return;
	
	out << "\n element group: " << fFEManager.ElementGroupNumber(this) + 1 << '\n';
	out << "\n element (in group): " << fElementCards.Position() + 1 << '\n';

	/* block data */
	const StringT& block_ID = ElementBlockID(fElementCards.Position());
	const ElementBlockDataT& block_data = BlockData(block_ID);

	/* model manager - block processor number */
	ModelManagerT* model = fFEManager.ModelManager();
	iArrayT elem_map(block_data.Dimension());
	model->ElementMap(block_ID, elem_map);

	/* block global number */
	const iArrayT* global_map = fFEManager.ElementMap(block_ID);

	/* report */
	out << "\n element (in partition block): " << elem_map[fElementCards.Position()] << '\n';
	if (global_map)
		out << " element (in global block): " << (*global_map)[fElementCards.Position()] << '\n';
	else
		out << " element (in global block): " << elem_map[fElementCards.Position()] << '\n';

	/* node number map */
	const iArrayT* node_map = fFEManager.NodeMap();
	iArrayT temp;

	out <<   "\n connectivity(x):\n";
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
  if (fConnectivities.Length() != fEqnos.Length())
    {
      cout << "ElementBaseT::SetElementCards length mismatch ";
      cout << "\n           element group: " << fFEManager.ElementGroupNumber(this) + 1;      
      cout << "\n fConnectivities length = " << fConnectivities.Length();
      cout << "\n          fEqnos length = " << fEqnos.Length() << endl;
      throw eSizeMismatch;
    }

	/* allocate */
	fElementCards.Allocate(fNumElements);

	/* loop over blocks to set pointers */
	int numberofnodes = fNodes->NumNodes();
	int count = 0;
	for (int i = 0; i < fBlockData.Length(); i++)
	{
		int dim = fBlockData[i].Dimension();
		int mat = fBlockData[i].MaterialID();
		const iArray2DT* blockconn = fConnectivities[i];
		iArray2DT& blockeqnos = fEqnos[i];

		if (blockconn->MajorDim() != blockeqnos.MajorDim())
		  {
		    cout << "ElementBaseT::SetElementCards length mismatch ";
		    cout << "\n   element group: " << fFEManager.ElementGroupNumber(this) + 1; 
		    cout << "\n           block: " << i+1;
		    cout << "\n  blockconn dim = " << blockconn->MajorDim() << " " << blockconn->MinorDim();
		    cout << "\n blockeqnos dim = " << blockeqnos.MajorDim() << " " << blockeqnos.MinorDim() << endl;
		    throw eSizeMismatch;
		  }

		for (int j = 0; j < dim; j++)
		{
			ElementCardT& element_card = fElementCards[count];
	
			/* material number */
			element_card.SetMaterialNumber(mat);

			/* set pointers */
			iArrayT& nodes = element_card.NodesX();
			blockconn->RowAlias(j, nodes);
			blockeqnos.RowAlias(j, element_card.Equations());
			
			/* check node numbers */
			int min, max;
			nodes.MinMax(min, max);
			if (min < 0 || max >= numberofnodes)
			{
				cout << "\n ElementBaseT::SetElementCards: nodes {" << min + 1
				     << "," << max + 1 << "} in element " << dim + 1 << "\n";
				cout <<   "     (" << j + 1 << " in block " <<  i + 1 << ") of group "
				     << fFEManager.ElementGroupNumber(this) + 1 << " are out of range" << endl;
				throw eBadInputValue;
			}

			count ++; /* next element */
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
