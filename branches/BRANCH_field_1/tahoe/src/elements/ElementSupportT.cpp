/* $Id: ElementSupportT.cpp,v 1.1.2.1 2002-04-27 01:32:26 paklein Exp $ */
#include "ElementSupportT.h"
#include "FEManagerT.h"
#include "NodeManagerPrimitive.h"

/* constructor */
inline ElementSupportT::ElementSupportT(void)
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

/* (re-)set the NodeManagerPrimitive */
void ElementSupportT::SetNodes(NodeManagerPrimitive* nodes)
{
	fNodes = nodes;
	if (nodes)
	{		
		fNumSD = nodes->NumSD();
	}
	else /* clear */
	{
		fNumSD = 0;
	}
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

int ElementSupportT::ElementGroupNumber(const ElementBaseT* element) const
{ 
	return FEManager().ElementGroupNumber(element); 
}

/* XDOF support */
XDOF_ManagerT& ElementSupportT::XDOF_Manager(void) const
{
	return Nodes();
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

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const
{
	FEManager().AssembleLHS(group, elMat, eqnos);
}

void ElementSupportT::AssembleRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const
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

/***********************************************************************
* Private
***********************************************************************/
