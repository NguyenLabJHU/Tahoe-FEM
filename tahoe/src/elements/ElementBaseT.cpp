/* $Id: ElementBaseT.cpp,v 1.46 2004-06-17 06:42:45 paklein Exp $ */
/* created: paklein (05/24/1996) */
#include "ElementBaseT.h"

#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "toolboxConstants.h"

#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "FieldT.h"
#include "eIntegratorT.h"
#endif

#include "LocalArrayT.h"

using namespace Tahoe;

/* constructor */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
ElementBaseT::ElementBaseT(const ElementSupportT& support, const FieldT& field):
	ParameterInterfaceT("element_base"),
	fSupport(support),
	fField(&field),
	fIntegrator(NULL),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{
	/* just cast it */
	fIntegrator = fSupport.eIntegrator(field);
}

ElementBaseT::ElementBaseT(const ElementSupportT& support):
	ParameterInterfaceT("element_base"),
	fSupport(support),
	fField(NULL),
	fIntegrator(NULL),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{

}
#else
ElementBaseT::ElementBaseT(ElementSupportT& support):
	ParameterInterfaceT("element_base"),
	fSupport(support),
	fIntegrator(NULL),
	fElementCards(0),
	fLHS(ElementMatrixT::kSymmetric)
{
	/* do nothing */
}
#endif

/* the index of this element group within the FEManagerT */
int ElementBaseT::ElementGroupNumber(void) const
{
	return ElementSupport().ElementGroupNumber(this);
}

/* destructor */
ElementBaseT::~ElementBaseT(void) {	}

/* allocates space and reads connectivity data */
void ElementBaseT::Initialize(void)
{
	/* set console variables */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int index = fSupport.ElementGroupNumber(this) + 1;
	StringT name;
	name.Append(index);
	name.Append("_element_group");
	iSetName(name);

	/* streams */
	ifstreamT& in = fSupport.Input();
	ostream&   out = fSupport.Output();

	/* control data */
	PrintControlData(out);

	/* element connectivity data */
	EchoConnectivityData(in, out);
#else
	EchoConnectivityData();
#endif

	/* dimension */
	int neq = NumElementNodes()*NumDOF();
	fLHS.Dimension(neq);	
	fRHS.Dimension(neq);
}

/* set the active elements */
void ElementBaseT::SetStatus(const ArrayT<StatusT>& status)
{
	if (status.Length() != NumElements()) ExceptionT::SizeMismatch();
	for (int i = 0; i < fElementCards.Length(); i++)
		fElementCards[i].Flag() = status[i];
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

#ifndef _FRACTURE_INTERFACE_LIBRARY_

/* the iteration number for the current time increment */
const int& ElementBaseT::IterationNumber(void) const
{
	return ElementSupport().IterationNumber(Group());
}

/* return true if the element contributes to the group */
bool ElementBaseT::InGroup(int group) const
{
	return Field().Group() == group;
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
void ElementBaseT::FormLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	try { LHSDriver(sys_type); }
	catch (ExceptionT::CodeT error)
	{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		cout << "\n ElementBaseT::FormLHS: " << fSupport.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fSupport.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{
			ostream& out = fSupport.Output();
		
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
#endif
		ExceptionT::Throw(error);
	}
}

void ElementBaseT::FormRHS(void)
{
	try { RHSDriver(); }
	catch (ExceptionT::CodeT error)
	{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		cout << "\n ElementBaseT::FormRHS: " << fSupport.Exception(error);
		cout << " in element " << fElementCards.Position() + 1 << " of group ";
		cout << fSupport.ElementGroupNumber(this) + 1 << ".\n";
		
		if (fElementCards.InRange())
		{	
			ostream& out = fSupport.Output();
	
			/* header */
			out << "\n ElementBaseT::FormRHS: caught exception " << error << '\n';
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
#endif
		ExceptionT::Throw(error);
	}
}

/* initialize/finalize time increment */
void ElementBaseT::InitStep(void) { }
void ElementBaseT::CloseStep(void) { }
GlobalT::RelaxCodeT ElementBaseT::ResetStep(void) 
{ 
	/* no relaxation */
	return GlobalT::kNoRelax;
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
	if (fConnectivities.Length() != fEqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* get local equations numbers */
		Field().SetLocalEqnos(*fConnectivities[i], fEqnos[i]);

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

#ifndef _FRACTURE_INTERFACE_LIBRARY_
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
#else
void ElementBaseT::ReadRestart(double* incomingData)
{
#pragma unused(incomingData)
	// Do nothing
}

void ElementBaseT::WriteRestart(double* outgoingData) const
{
#pragma unused(outgoingData)
	// Do nothing
}
#endif

/* returns 1 if DOF's are interpolants of the nodal values */
int ElementBaseT::InterpolantDOFs(void) const { return 1; }

void ElementBaseT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
#if __option(extended_errorcheck)
	if (nodes.Length() != DOFs.MajorDim() ||
	    DOFs.MinorDim() != NumDOF()) ExceptionT::SizeMismatch("ElementBaseT::NodalDOFs");
#endif

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	const dArray2DT& all_DOFs = Field()[0]; // displacements
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
	const char caller[] = "ElementBaseT::ElementBlockID";
	if (element < 0 || element >= NumElements())
		ExceptionT::OutOfRange(caller, "element number %d is out of range {0,%d}", element, NumElements() - 1);
	
	bool found = false;
	for (int i = 0; i < fBlockData.Length(); i++)
		if (element >= fBlockData[i].StartNumber() &&
		    element <  fBlockData[i].StartNumber() + fBlockData[i].Dimension())
			return fBlockData[i].ID();

	if (!found)
		ExceptionT::GeneralFail(caller, "could not resolve block ID for element %d", element);

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
		const int* p = elemnodes.Pointer();
		for (int n=0; n < elemnodes.Length(); n++)
			if (weight[*p] < base_weight) 
				weight[*p] = base_weight;
	}
}

/* array of nodes used by the element group */
void ElementBaseT::NodesUsed(ArrayT<int>& nodes_used) const
{
	int num_blocks = fBlockData.Length();
	if (num_blocks == 0)
		ExceptionT::GeneralFail("ElementBaseT::NodesUsed", "fBlockData is not set");

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
		const int *pc = conn->Pointer();
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

/* add the element group's contribution to the lumped (scalar) mass of the given nodes */
void ElementBaseT::LumpedMass(const iArrayT& nodes, dArrayT& mass) const
{
	if (nodes.Length() != mass.Length())
		ExceptionT::SizeMismatch("ElementBaseT::LumpedMass");
}

/* contribution to the nodal residual forces */
const dArray2DT& ElementBaseT::InternalForce(int group)
{
#pragma unused(group)
	ExceptionT::GeneralFail("ElementBaseT::InternalForce", "not implemented");
	return ElementSupport().CurrentCoordinates(); /* dummy */
}

/* describe the parameters needed by the interface */
void ElementBaseT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* associated fields */
	list.AddParameter(ParameterT::String, "field_name");
}

/*information about subordinate parameter lists */
void ElementBaseT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	sub_list.AddSub("element_block", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ElementBaseT::NewSub(const StringT& list_name) const
{
	if (list_name == "element_block")
		return new ElementBlockDataT;
	else
		return ParameterInterfaceT::NewSub(list_name);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* map the element numbers from block to group numbering */
void ElementBaseT::BlockToGroupElementNumbers(iArrayT& elems, const StringT& block_ID) const
{
	const char caller[] = "ElementBaseT::BlockToGroupElementNumbers";

	/* block information */
	const ElementBlockDataT& block_data = BlockData(block_ID);

	/* check */
	int min, max;
	elems.MinMax (min, max);
	if (min < 0 || max > block_data.Dimension())
		ExceptionT::BadInputValue(caller, "element numbers {%d,%d} are out of range",
			min, max);
	
	/* map to group numbering is just a shift */
	elems += block_data.StartNumber();
}

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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	fSupport.AssembleRHS(fField->Group(), fRHS, CurrentElement().Equations());
#else
	fSupport.AssembleRHS(fElementCards.Position(),fRHS,CurrentElement().Equations());
#endif
}

void ElementBaseT::AssembleLHS(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	fSupport.AssembleLHS(fField->Group(), fLHS, CurrentElement().Equations());
#else	
	fSupport.AssembleLHS(fElementCards.Position(),fLHS,CurrentElement().Equations());
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
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

#else

void ElementBaseT::EchoConnectivityData(void)
{	
	/* read */
	ReadConnectivity();

#endif

	/* derived dimensions */
	int neq = NumElementNodes()*NumDOF();
	fEqnos.Dimension(fBlockData.Length());
	for (int be=0; be < fEqnos.Length(); be++)
	  {
	    int numblockelems = fConnectivities[be]->MajorDim();
	    fEqnos[be].Dimension(numblockelems, neq);
	    fEqnos[be] = -1;
	  }

	/* set pointers in element cards */
	SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* write */
	WriteConnectivity(out);
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* resolve input format types */
void ElementBaseT::ReadConnectivity(ifstreamT& in, ostream& out)
{
#pragma unused(out)
#else
void ElementBaseT::ReadConnectivity(void)
{
#endif

	/* read from parameter file */
	ArrayT<StringT> elem_ID;
	iArrayT matnums;
	ModelManagerT& model = fSupport.Model();
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	model.ElementBlockList(in, elem_ID, matnums);
#else
	/* For Sierra, can't use input stream */
	elem_ID.Dimension(1);
	matnums.Dimension(1);
	elem_ID[0] = fSupport.BlockID();
	/*Might have to generalize this later*/
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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	    model.ElementGroupDimensions(elem_ID[b], num_elems, num_nodes);
#else
		num_elems = fSupport.NumElements();
		num_nodes = model.NumNodes(); 
#endif

	    /* set if unset */
	    if (nen == 0) nen = num_nodes;
	    
	    /* consistency check */
	    if (num_nodes != 0 && nen != num_nodes)
			ExceptionT::BadInputValue("ElementBaseT::ReadConnectivity",
				"minor dimension %d of block %d does not match previous %d", num_nodes, b+1, nen);
	    
	    /* store block data */
	    fBlockData[b].Set(elem_ID[b], elem_count, num_elems, matnums[b] - 1); // offset

	    /* increment element count */
	    elem_count += num_elems;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
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

#ifndef _FRACTURE_INTERFACE_LIBRARY_
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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
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
#endif
		ExceptionT::BadInputValue("ElementBaseT::BlockData");
	}

	/* return */
	return fBlockData[block_num];
}

/* write all current element information to the stream */
void ElementBaseT::CurrElementInfo(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
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
	const ArrayT<int>* node_map = fSupport.NodeMap();
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
#else
#pragma unused(out)
#endif 
}

/* set element cards array */
void ElementBaseT::SetElementCards(
	const ArrayT<ElementBlockDataT>& block_data, 
	const ArrayT<const iArray2DT*>& connectivities,		
	const ArrayT<iArray2DT>& eqnos, 
	AutoArrayT<ElementCardT>& element_cards) const
{
	const char caller[] = "ElementBaseT::SetElementCards";
	if (connectivities.Length() != eqnos.Length())
    {
#ifndef _FRACTURE_INTERFACE_LIBRARY_
      cout << "ElementBaseT::SetElementCards length mismatch ";
      cout << "\n           element group: " << fSupport.ElementGroupNumber(this) + 1;      
      cout << "\n connectivities length = " << connectivities.Length();
      cout << "\n         eqnos length = " << eqnos.Length() << endl;
#endif
		ExceptionT::SizeMismatch(caller);
    }

	/* loop over blocks to set pointers */
	int numberofnodes = fSupport.NumNodes();
	int count = 0;
	for (int i = 0; i < block_data.Length(); i++)
	{
		int dim = block_data[i].Dimension();
		int mat = block_data[i].MaterialID();
		const iArray2DT* blockconn = connectivities[i];
		const iArray2DT& blockeqnos = eqnos[i];

		if (blockconn->MajorDim() != blockeqnos.MajorDim())
		  {
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		    cout << "ElementBaseT::SetElementCards length mismatch ";
		    cout << "\n   element group: " << fSupport.ElementGroupNumber(this) + 1; 
		    cout << "\n           block: " << i+1;
		    cout << "\n  blockconn dim = " << blockconn->MajorDim() << " " << blockconn->MinorDim();
		    cout << "\n blockeqnos dim = " << blockeqnos.MajorDim() << " " << blockeqnos.MinorDim() << endl;
#endif
		    ExceptionT::SizeMismatch(caller);
		  }

		for (int j = 0; j < dim; j++)
		{
			ElementCardT& element_card = element_cards[count];
	
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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
				cout << "\n ElementBaseT::SetElementCards: nodes {" << min + 1
				     << "," << max + 1 << "} in element " << dim + 1 << "\n";
				cout <<   "     (" << j + 1 << " in block " <<  i + 1 << ") of group "
				     << fSupport.ElementGroupNumber(this) + 1 << " are out of range" << endl;
#endif
				ExceptionT::BadInputValue(caller);
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
