/* $Id: ElementBaseT.cpp,v 1.24.2.1 2002-11-13 08:37:52 paklein Exp $ */
/* created: paklein (05/24/1996) */

#include "ElementBaseT.h"

#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "ModelManagerT.h"
#include "fstreamT.h"
#include "toolboxConstants.h"

#ifndef _SIERRA_TEST_
#include "FieldT.h"
#include "eControllerT.h"
#endif

#include "LocalArrayT.h"

using namespace Tahoe;

/* constructor */
#ifndef _SIERRA_TEST_
ElementBaseT::ElementBaseT(const ElementSupportT& support, const FieldT& field):
	fSupport(support),
	fField(field),
	fController(NULL),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{
	/* just cast it */
	fController = fSupport.eController(field);
}
#else
ElementBaseT::ElementBaseT(const ElementSupportT& support):
	fSupport(support),
	fController(NULL),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{
	/* do nothing */
}
#endif

/* destructor */
ElementBaseT::~ElementBaseT(void) {	}

/* allocates space and reads connectivity data */
void ElementBaseT::Initialize(void)
{
	/* set console variables */
	int index = fSupport.ElementGroupNumber(this) + 1;
	StringT name;
	name.Append(index);
	name.Append("_element_group");
	iSetName(name);

	/* streams */
	ifstreamT& in = fSupport.Input();
	ostream&   out = fSupport.Output();
#ifndef _SIERRA_TEST_
	/* control data */
	PrintControlData(out);
#endif
	/* element connectivity data */
	EchoConnectivityData(in, out);

	/* dimension */
	int neq = NumElementNodes()*NumDOF();
	fLHS.Dimension(neq);	
	fRHS.Dimension(neq);
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

/* form of tangent matrix - symmetric by default */
GlobalT::SystemTypeT ElementBaseT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

#ifndef _SIERRA_TEST_
/* the iteration number for the current time increment */
const int& ElementBaseT::IterationNumber(void) const
{
	return ElementSupport().IterationNumber(Group());
}
#endif

/* collect the list of element block ID's used by the element group */
void ElementBaseT::ElementBlockIDs(ArrayT<StringT>& IDs) const
{
	IDs.Dimension(fBlockData.Length());
	for (int i = 0; i < IDs.Length(); i++)
		IDs[i] = fBlockData[i].ID();
}

/* solution calls */
void ElementBaseT::FormLHS(void)
{
	try { LHSDriver(); }
	catch (ExceptionT::CodeT error)
	{
		cout << "\n ElementBaseT::FormLHS: " << fSupport.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fSupport.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{
#ifndef _SIERRA_TEST_		
			ostream& out = fSupport.Output();
#else
			ostream& out = cout;
#endif
		
			/* header */
			out << "\n ElementBaseT::FormLHS: caught exception " << error << '\n';
			out <<   "      Time: " << fSupport.Time() << '\n';
			out <<   "      Step: " << fSupport.StepNumber() << '\n';
			out <<   " Time step: " << fSupport.TimeStep() << '\n';
		
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
	catch (ExceptionT::CodeT error)
	{
		cout << "\n ElementBaseT::FormRHS: " << fSupport.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fSupport.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{
#ifndef _SIERRA_TEST_		
			ostream& out = fSupport.Output();
	
			/* header */
			out << "\n ElementBaseT::FormRHS: caught exception " << error << '\n';
			out <<   "      Time: " << fSupport.Time() << '\n';
			out <<   "      Step: " << fSupport.StepNumber() << '\n';
			out <<   " Time step: " << fSupport.TimeStep() << '\n';
		
			/* write current element information to main out */
			CurrElementInfo(out);
			cout << "     See output file for current element information\n";
#endif
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
void ElementBaseT::ResetStep(void) { }

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
	if (fConnectivities.Length() != fEqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

#ifndef _SIERRA_TEST_
	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* get local equations numbers */
		fField.SetLocalEqnos(*fConnectivities[i], fEqnos[i]);

		/* add to list of equation numbers */
		eq_1.Append(&fEqnos[i]);
	}
#else
#pragma unused(eq_1)
#endif
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

void ElementBaseT::ReadRestart(istream& in)
{
	/* stream check */
	if (!in.good()) throw ExceptionT::kGeneralFail;
}

void ElementBaseT::WriteRestart(ostream& out) const
{
	/* stream check */
	if (!out.good()) throw ExceptionT::kGeneralFail;
}

/* returns 1 if DOF's are interpolants of the nodal values */
int ElementBaseT::InterpolantDOFs(void) const { return 1; }

void ElementBaseT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
#if __option(extended_errorcheck)
	if (nodes.Length() != DOFs.MajorDim() ||
	    DOFs.MinorDim() != NumDOF()) throw ExceptionT::kSizeMismatch;

#endif

#ifndef _SIERRA_TEST_
	const dArray2DT& all_DOFs = fField[0]; // displacements
	DOFs.RowCollect(nodes, all_DOFs); // no check of nodes used

//NOTE - This function is added only for completeness. If the
//       DOF's are interpolant, there should be no reason to
//       collect the displacements in this way. Use SetLocalU
#else
#pragma unused(nodes)
#pragma unused(DOFs)
#endif
}

/* block ID for the specified element */
const StringT& ElementBaseT::ElementBlockID(int element) const
{
	if (element < 0 || element >= NumElements()) {
		cout << "\n ElementBaseT::ElementBlockID: element number " << element << " is out of range {0,"
		    << NumElements() - 1 << "}" << endl;
		throw ExceptionT::kOutOfRange;
	}
	
	bool found = false;
	for (int i = 0; i < fBlockData.Length(); i++)
		if (element >= fBlockData[i].StartNumber() &&
		    element <  fBlockData[i].StartNumber() + fBlockData[i].Dimension())
			return fBlockData[i].ID();

	if (!found) {
		cout << "\n ElementBaseT::ElementBlockID: could not resolve block ID for element "
		     << element << endl;
		throw ExceptionT::kGeneralFail;
	}
	return fBlockData[0].ID(); /* dummy */
}

/* weight the computational effort of every node */
void ElementBaseT::WeightNodalCost(iArrayT& weight) const
{
	int base_weight = 1;
	int nel = NumElements();
	for (int i=0; i < nel; i++)
	{
		const iArrayT& elemnodes = fElementCards[i].NodesX();
		int* p = elemnodes.Pointer();
		for (int n=0; n < elemnodes.Length(); n++)
			if (weight[*p] < base_weight) 
				weight[*p] = base_weight;
	}
}

/* moved from protected to public by HSP on 7-19-02 */
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
	nodes_used.Dimension(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodes_used[dex++] = j + min;
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
#ifndef _SIERRA_TEST_
	fSupport.AssembleRHS(fField.Group(), fRHS, CurrentElement().Equations());
#else
	fSupport.AssembleRHS(fElementCards.Position(),fRHS,CurrentElement().Equations());
#endif
}

void ElementBaseT::AssembleLHS(void) const
{
#ifndef _SIERRA_TEST_
	fSupport.AssembleLHS(fField.Group(), fLHS, CurrentElement().Equations());
#endif
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
#ifndef _SIERRA_TEST_	
	out << "\n Element Connectivity:\n";
#endif	
	/* read */
	ReadConnectivity(in, out);

	/* derived dimensions */
	int neq = NumElementNodes()*NumDOF();
	fEqnos.Dimension(fBlockData.Length());
	for (int be=0; be < fEqnos.Length(); be++)
	  {
	    int numblockelems = fConnectivities[be]->MajorDim();
	    fEqnos[be].Dimension(numblockelems, neq);
	  }

	/* set pointers in element cards */
	SetElementCards();

#ifndef _SIERRA_TEST_
	/* write */
	WriteConnectivity(out);
#endif
}

/* resolve input format types */
void ElementBaseT::ReadConnectivity(ifstreamT& in, ostream& out)
{
#pragma unused(out)
#ifdef _SIERRA_TEST_
#pragma unused(in)
#endif

	/* read from parameter file */
	ArrayT<StringT> elem_ID;
	iArrayT matnums;
	ModelManagerT& model = fSupport.Model();
#ifndef _SIERRA_TEST_
	model.ElementBlockList(in, elem_ID, matnums);
#else
	/* For Sierra, can't use input stream */
	elem_ID.Dimension(1);
	matnums.Dimension(1);
	elem_ID = "1";
	matnums = 1; 
#endif

	/* allocate block map */
	int num_blocks = elem_ID.Length();
	fBlockData.Dimension(num_blocks);
	fConnectivities.Allocate (num_blocks);

	/* read from parameter file */
	int elem_count = 0;
	int nen = 0;
	for (int b=0; b < num_blocks; b++)
	{
	    /* check number of nodes */
	    int num_elems, num_nodes;
#ifndef _SIERRA_TEST_
	    model.ElementGroupDimensions(elem_ID[b], num_elems, num_nodes);
#else
		num_elems = fSupport.NumElements();
		num_nodes = model.NumNodes(); 
#endif

	    /* set if unset */
	    if (nen == 0) nen = num_nodes;
	    
	    /* consistency check */
	    if (num_nodes != 0 && nen != num_nodes)
		{
#ifndef _SIERRA_TEST_
			cout << "\n ElementBaseT::ReadConnectivity: minor dimension "
                 << num_nodes << " of block " << b+1 << '\n';
			cout <<   "     does not match dimension of previous blocks "
                 << nen << endl;
#endif                 
			throw ExceptionT::kBadInputValue;
		}
	    
	    /* store block data */
	    fBlockData[b].Set(elem_ID[b], elem_count, num_elems, matnums[b] - 1); // offset

	    /* increment element count */
	    elem_count += num_elems;

#ifndef _SIERRA_TEST_
	    /* load connectivity from database into model manager */
	    model.ReadConnectivity(elem_ID[b]);
#endif

	    /* set pointer to connectivity list */
	    fConnectivities[b] = model.ElementGroupPointer(elem_ID[b]);
	}

	/* connectivities came back empty */
	if (nen == 0) nen = DefaultNumElemNodes();
	for (int i = 0; i < fConnectivities.Length(); i++)
		if (fConnectivities[i]->MinorDim() == 0)
		{
			/* not really violating const-ness */
			iArray2DT* connects = const_cast<iArray2DT*>(fConnectivities[i]);
			connects->Dimension(0, nen);
		}
	  
	/* set dimensions */
	fElementCards.Dimension(elem_count);
}

/* resolve output formats */
void ElementBaseT::WriteConnectivity(ostream& out) const
{	
	out << " Number of elements. . . . . . . . . . . . . . . = " << NumElements() << '\n';

	/* write dimensions of blocks */
	out << " Block dimensions:\n";
	out << setw(kIntWidth) << "ID"
	    << setw(kIntWidth) << "size" << '\n';
	for (int i = 0; i < fBlockData.Length(); i++)
		out << setw(kIntWidth) << fBlockData[i].ID()
		    << setw(kIntWidth) << fBlockData[i].Dimension() << '\n';
	out << endl;

	/* verbose output */
	if (fSupport.PrintInput())
	{
		/* write header */
		out << setw(kIntWidth) << "no.";
		out << setw(kIntWidth) << "mat.";
		int nen = NumElementNodes();
		for (int j = 1; j <= ((nen < 9) ? nen : 8); j++)
		{
			int numwidth = (j < 10) ? 1 : ((j < 100) ? 2 : 3);		
			out << setw(kIntWidth - (numwidth + 1)) << "n[";
			out << j << "]";
		}
		out << endl;
				
		/* write material number and connectivity */
		iArrayT nodesX(nen);
		for (int i = 0; i < NumElements(); i++)
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

#if 0
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
	localconnects.Dimension(NumElements(), NumElementNodes());
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
#endif

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
		cout <<   "     element group " << fSupport.ElementGroupNumber(this) + 1;
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
		throw ExceptionT::kBadInputValue;
	}

	/* return */
	return fBlockData[block_num];
}

/* write all current element information to the stream */
void ElementBaseT::CurrElementInfo(ostream& out) const
{
#pragma unused(out)
#ifndef _SIERRA_TEST_
	if (!fElementCards.InRange()) return;
	
	out << "\n element group: " << fSupport.ElementGroupNumber(this) + 1 << '\n';
	out << "\n element (in group): " << fElementCards.Position() + 1 << '\n';

	/* block data */
	const StringT& block_ID = ElementBlockID(fElementCards.Position());
	const ElementBlockDataT& block_data = BlockData(block_ID);
	int local_el_number = fElementCards.Position() - block_data.StartNumber();

	/* model manager - block processor number */
	ModelManagerT& model = fSupport.Model();
	iArrayT elem_map(block_data.Dimension());
	model.ElementMap(block_ID, elem_map);

	/* block global number */
	const iArrayT* global_map = fSupport.ElementMap(block_ID);

	/* report */
	out << "\n element (in partition block): " << elem_map[local_el_number] << '\n';
	if (global_map)
		out << " element (in global block): " << (*global_map)[local_el_number] << '\n';
	else
		out << " element (in global block): " << elem_map[local_el_number] << '\n';

	/* node number map */
	const iArrayT* node_map = fSupport.NodeMap();
	iArrayT temp;

	out <<   "\n connectivity(x):\n";
	if (node_map)
	{
		const iArrayT& nodes_X = CurrentElement().NodesX();
		temp.Dimension(nodes_X.Length());
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
		temp.Dimension(nodes_U.Length());
		for (int i = 0; i < nodes_U.Length(); i++)
			temp[i] = (*node_map)[nodes_U[i]];
	}
	else
		temp = CurrentElement().NodesU();
	temp++;
	out << temp.wrap(4) << '\n';

	out <<   " equations:\n";
	out << (CurrentElement().Equations()).wrap(4) << '\n';
#endif // ndef _SIERRA_TEST_
}

/* set element cards array */
void ElementBaseT::SetElementCards(void)
{
  if (fConnectivities.Length() != fEqnos.Length())
    {
      cout << "ElementBaseT::SetElementCards length mismatch ";
      cout << "\n           element group: " << fSupport.ElementGroupNumber(this) + 1;      
      cout << "\n fConnectivities length = " << fConnectivities.Length();
      cout << "\n          fEqnos length = " << fEqnos.Length() << endl;
      throw ExceptionT::kSizeMismatch;
    }

	/* allocate */
	//fElementCards.Dimension(fNumElements);

	/* loop over blocks to set pointers */
	int numberofnodes = fSupport.NumNodes();
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
		    cout << "\n   element group: " << fSupport.ElementGroupNumber(this) + 1; 
		    cout << "\n           block: " << i+1;
		    cout << "\n  blockconn dim = " << blockconn->MajorDim() << " " << blockconn->MinorDim();
		    cout << "\n blockeqnos dim = " << blockeqnos.MajorDim() << " " << blockeqnos.MinorDim() << endl;
		    throw ExceptionT::kSizeMismatch;
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
				     << fSupport.ElementGroupNumber(this) + 1 << " are out of range" << endl;
				throw ExceptionT::kBadInputValue;
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