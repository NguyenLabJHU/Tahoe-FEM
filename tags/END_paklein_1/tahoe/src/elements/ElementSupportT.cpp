/* $Id: ElementSupportT.cpp,v 1.4.6.1 2002-10-17 04:28:48 paklein Exp $ */
#include "ElementSupportT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "eControllerT.h"
#include "nControllerT.h"
#include "FieldT.h"

/* constructor */

using namespace Tahoe;

ElementSupportT::ElementSupportT(void)
{
	/* clear */
	SetFEManager(NULL);
}

/* (re-)set the FEManagerT */
void ElementSupportT::SetFEManager(FEManagerT* fe)
{
	fFEManager = fe;
	if (fe)
	{
		fAnalysis = fe->Analysis();
		fRunState = &(fe->RunState());

		/* set nodal information */
		SetNodes(fe->NodeManager());
	}
	else
	{
		fAnalysis = GlobalT::kNoAnalysis;
		fRunState = NULL;

		/* clear nodal information */
		SetNodes(NULL);
	}
}

/* (re-)set the NodeManagerT */
void ElementSupportT::SetNodes(NodeManagerT* nodes)
{
	fNodes = nodes;
	if (nodes)
	{		
		fNumSD = nodes->NumSD();
		fNumNodes = nodes->NumNodes();
	}
	else /* clear */
	{
		fNumSD = 0;
		fNumNodes = 0;
	}
}

/* Tahoe version string */
const StringT& ElementSupportT::Version(void) const
{
	return FEManager().Version();
}

bool ElementSupportT::PrintInput(void) const
{
	return FEManager().PrintInput();
}

const dArray2DT& ElementSupportT::InitialCoordinates(void) const
{
	return Nodes().InitialCoordinates();
}

const dArray2DT& ElementSupportT::CurrentCoordinates(void) const
{
	return Nodes().CurrentCoordinates();
}

void ElementSupportT::RegisterCoordinates(LocalArrayT& array) const
{
	Nodes().RegisterCoordinates(array);
}

/* return a  schedule function */
const ScheduleT* ElementSupportT::Schedule(int num) const
{
	return FEManager().Schedule(num);
}

const int& ElementSupportT::IterationNumber(int group) const 
{ 
	return FEManager().IterationNumber(group); 
}

const char* ElementSupportT::Exception(ExceptionT::CodeT exception) const
{
	return ExceptionT::ToString(exception);
}

int ElementSupportT::ElementGroupNumber(const ElementBaseT* element) const
{ 
	return FEManager().ElementGroupNumber(element); 
}

const double& ElementSupportT::Time(void) const
{
	return FEManager().Time();
}

const double& ElementSupportT::TimeStep(void) const
{
	return FEManager().TimeStep();
}

const int& ElementSupportT::StepNumber(void) const
{
	return FEManager().StepNumber();
}

const int& ElementSupportT::NumberOfSteps(void) const
{
	return FEManager().NumberOfSteps();
}

/* the element group at the specified index in the element list */
ElementBaseT& ElementSupportT::ElementGroup(int index) const
{
	ElementBaseT* element = FEManager().ElementGroup(index);
	if (!element) throw ExceptionT::kGeneralFail;
	return *element;
}

/* geometry information */
ModelManagerT& ElementSupportT::Model(void) const
{
	ModelManagerT* model = FEManager().ModelManager();
	if (!model) throw ExceptionT::kGeneralFail;
	return *model;
}

/* XDOF support */
XDOF_ManagerT& ElementSupportT::XDOF_Manager(void) const
{
	return Nodes();
}

/* node number map. returns NULL if there is not map */
const iArrayT* ElementSupportT::NodeMap(void) const
{
	return FEManager().NodeMap();
}

/* return a pointer to the field */
const FieldT* ElementSupportT::Field(const char* name) const
{
	return Nodes().Field(name);
}

/* return the element controller appropriate for the given field */
const eControllerT* ElementSupportT::eController(const FieldT& field) const
{
	const nControllerT& n_cont = field.nController();
	const eControllerT* e_cont = dynamic_cast<const eControllerT*>(&n_cont);
	return e_cont;
}
	
/* element number map for the given block ID */
const iArrayT* ElementSupportT::ElementMap(const StringT& block_ID) const
{
	return FEManager().ElementMap(block_ID);
}

/* MP */
int ElementSupportT::Size(void) const { return FEManager().Size(); }
int ElementSupportT::Rank(void) const { return FEManager().Rank(); }
void ElementSupportT::IncomingNodes(iArrayT& nodes_in) const
{
	FEManager().IncomingNodes(nodes_in);
}
void ElementSupportT::OutgoingNodes(iArrayT& nodes_out) const
{
	FEManager().OutgoingNodes(nodes_out);
}
void ElementSupportT::SendExternalData(const dArray2DT& all_out_data) const
{
	FEManager().SendExternalData(all_out_data);
}
void ElementSupportT::RecvExternalData(dArray2DT& external_data) const
{
	FEManager().RecvExternalData(external_data);
}

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& eqnos) const
{
	FEManager().AssembleLHS(group, elMat, eqnos);
}

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos) const
{
	FEManager().AssembleLHS(group, elMat, row_eqnos, col_eqnos);
}

void ElementSupportT::AssembleRHS(int group, const dArrayT& elRes, 
	const nArrayT<int>& eqnos) const
{
	FEManager().AssembleRHS(group, elRes, eqnos);
}

/* initialize work space to the number of values to be averaged */
void ElementSupportT::ResetAverage(int n_values) const
{
	Nodes().ResetAverage(n_values);
}

/* assemble values */
void ElementSupportT::AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const
{
	Nodes().AssembleAverage(nodes, vals);
}

/* average assembled values */
const dArray2DT& ElementSupportT::OutputAverage(void) const
{
	return Nodes().OutputAverage();
}

/* return averaged values for the nodes with assembled values */
void ElementSupportT::OutputUsedAverage(dArray2DT& average_values) const
{
	Nodes().OutputUsedAverage(average_values);
}

ifstreamT& ElementSupportT::Input(void) const
{
	return FEManager().Input();
}

ofstreamT& ElementSupportT::Output(void) const
{
	return FEManager().Output();
}

int ElementSupportT::RegisterOutput(const OutputSetT& output_set) const
{
	return FEManager().RegisterOutput(output_set);
}

void ElementSupportT::WriteOutput(int ID, const dArray2DT& n_values, 
	const dArray2DT& e_values) const
{
	FEManager().WriteOutput(ID, n_values, e_values);
}

/***********************************************************************
* Private
***********************************************************************/
