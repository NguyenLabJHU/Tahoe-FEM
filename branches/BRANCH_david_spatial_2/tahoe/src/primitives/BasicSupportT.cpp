/* $Id: BasicSupportT.cpp,v 1.6.8.2 2005-06-09 03:28:47 paklein Exp $ */
#include "BasicSupportT.h"

#include "dArray2DT.h"

#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "FEManagerT.h"
#include "TimeManagerT.h"
#include "CommManagerT.h"
#include "NodeManagerT.h"
#include "eIntegratorT.h"
#include "nIntegratorT.h"
#include "FieldT.h"
#include "CommunicatorT.h"
#else
#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dMatrixT.h"
#include "ElementMatrixT.h"
#include "IOBaseT.h"
#endif

using namespace Tahoe;

/* constructor */
BasicSupportT::BasicSupportT(void):
	fFEManager(NULL),
	fNodeManager(NULL),
	fTimeManager(NULL),
	fModelManager(NULL),
	fCommManager(NULL),
	fCommunicator(NULL),
	fRunState(NULL),
	fNumSD(-1)
{
	/* clear */
	SetFEManager(NULL);

#ifdef _FRACTURE_INTERFACE_LIBRARY_
	fNumSD = 3;
	fTimeStep = 0.0;
	fIterationNumber = 0;
	fCurrentCoordinates = NULL;
	fInitialCoordinates = NULL;
	fRHS = NULL;
	fLHS = NULL;
#endif
}

/* (re-)set the FEManagerT */
void BasicSupportT::SetFEManager(const FEManagerT* fe)
{
	fFEManager = fe;
	if (!fe)
	{
		fRunState = NULL;

		/* clear nodal information */
		SetNodeManager(NULL);

		/* clear pointers */
		fModelManager = NULL;
		fTimeManager = NULL;
		fCommManager = NULL;
		fCommunicator = NULL;
	}
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	else
	{
		fRunState = &(fFEManager->RunState());

		/* set nodal information */
		SetNodeManager(fFEManager->NodeManager());

		/* set model manager */
		fModelManager = fFEManager->ModelManager();

		/* set time manager */
		fTimeManager = fFEManager->TimeManager();

		/* set comm manager */
		fCommManager = fFEManager->CommManager();
		
		/* set communicator */
		if (fCommManager)
			fCommunicator = &(fCommManager->Communicator());
		else
			fCommunicator = NULL;
	}
#endif
}

/* (re-)set the NodeManagerT */
void BasicSupportT::SetNodeManager(NodeManagerT* nodes) { fNodeManager = nodes; }

void BasicSupportT::SetNumSD(int nsd)
{
	/* error check */
	if (nsd < -1 || nsd == 0)
		ExceptionT::GeneralFail("BasicSupportT::SetNumSD", "expecting n > 0 or n = -1, not n = %d", nsd);

	/* set cached value */
	fNumSD = nsd;
}

/* Tahoe version string */
const char* BasicSupportT::Version(void) const 
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return "0.1";
#else
	return FEManager().Version();
#endif
}

bool BasicSupportT::PrintInput(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return false;
#else
	return FEManager().PrintInput();
#endif
}

GlobalT::LoggingT BasicSupportT::Logging(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return GlobalT::kSilent;
#else
	return FEManager().Logging();
#endif
}

const dArray2DT& BasicSupportT::InitialCoordinates(void) const {
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	if (!fInitialCoordinates) ExceptionT::GeneralFail("BasicSupportT::InitialCoordinates", "pointer not set");
	return *fInitialCoordinates;
#else
	return NodeManager().InitialCoordinates();
#endif
}

const dArray2DT& BasicSupportT::CurrentCoordinates(void) const {
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	if (!fCurrentCoordinates) ExceptionT::GeneralFail("BasicSupportT::CurrentCoordinates", "pointer not set");
	return *fCurrentCoordinates;
#else
	return NodeManager().CurrentCoordinates();
#endif
}

void BasicSupportT::RegisterCoordinates(LocalArrayT& array) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
    switch (array.Type())
    {
    	case LocalArrayT::kInitCoords: 
    	{
    		array.SetGlobal(InitialCoordinates());
    		break;
    	}
    	case LocalArrayT::kCurrCoords:
    	{
    		array.SetGlobal(CurrentCoordinates());
    		break;
    	}
    	default:
			ExceptionT::GeneralFail("BasicSupportT::RegisterCoordinates", "unsupported type %d", array.Type());
     }
#else
	NodeManager().RegisterCoordinates(array);
#endif
}

/* return a  schedule function */
const ScheduleT* BasicSupportT::Schedule(int num) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#pragma unused(num)
	return NULL;
#else
	return FEManager().Schedule(num);
#endif
}

/* return the iteration number for the current solver group */
int BasicSupportT::IterationNumber(void) const
{ 
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fIterationNumber;
#else
	return FEManager().IterationNumber();
#endif
}

const int& BasicSupportT::IterationNumber(int group) const 
{ 
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#pragma unused(group)
	return fIterationNumber;
#else
	return FEManager().IterationNumber(group);
#endif
}

/* the group number being solved or -1 if not defined */
int BasicSupportT::CurrentGroup(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return -1;
#else
	return FEManager().CurrentGroup();
#endif
}

const double& BasicSupportT::Time(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fTimeStep;
#else
	return TimeManager().Time();
#endif
}

const double& BasicSupportT::TimeStep(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fTimeStep;
#else	
	return TimeManager().TimeStep();
#endif
}

const int& BasicSupportT::StepNumber(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fIterationNumber;
#else
	return TimeManager().StepNumber();
#endif
}

const int& BasicSupportT::NumberOfSteps(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return fIterationNumber;
#else
	return TimeManager().NumberOfSteps();
#endif
}

/* return a pointer to the field */
const FieldT* BasicSupportT::Field(const char* name) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#pragma unused(name)
	return NULL;
#else
	return NodeManager().Field(name);
#endif
}

/* node number map. returns NULL if there is not a map */
const ArrayT<int>* BasicSupportT::NodeMap(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return NULL;
#else
	return FEManager().NodeMap();	
#endif
}

/* element number map for the given block ID */
const iArrayT* BasicSupportT::ElementMap(const StringT& block_ID) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#pragma unused(block_ID)
	return NULL;	
#else
	return FEManager().ElementMap(block_ID);
#endif
}

int BasicSupportT::Size(void) const 
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return 1;
#else
	return Communicator().Size(); 	
#endif
}

int BasicSupportT::Rank(void) const 
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return 0;
#else
	return Communicator().Rank();	
#endif 
}

/* the local node to home processor map */
const ArrayT<int>* BasicSupportT::ProcessorMap(void) const 
{ 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	if (fFEManager)
		return fFEManager->ProcessorMap(); 
	else
#endif
	return NULL;
}

const ArrayT<int>* BasicSupportT::ExternalNodes(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	if (fCommManager) 
		return fCommManager->ExternalNodes();
	else
#endif
	return NULL;
}

/* list of nodes owned by this processor or NULL if \e all nodes are owned */
const ArrayT<int>* BasicSupportT::PartitionNodes(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	if (fCommManager) 
		return fCommManager->PartitionNodes();
	else
#endif
	return NULL;
}

void BasicSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleLHS(group, elMat, eqnos);
#else
#pragma unused(group)
#pragma unused(eqnos)
#pragma unused(elMat)
#endif
}

void BasicSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleLHS(group, elMat, row_eqnos, col_eqnos);
#else
#pragma unused(group)
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
#endif
}

void BasicSupportT::AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, 
	const nArrayT<int>& eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleLHS(group, diagonal_elMat, eqnos);
#else
#pragma unused(group)
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
#endif
}

void BasicSupportT::AssembleRHS(int group, const nArrayT<double>& elRes, 
	const nArrayT<int>& eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleRHS(group, elRes, eqnos);
#else
#pragma unused(group)
#pragma unused(eqnos)
#pragma unused(elRes)
#pragma unused(eqnos)
#endif
}

/* the residual for the given group */
const dArrayT& BasicSupportT::RHS(int group) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().RHS(group);
#else
#pragma unused(group)
	return *fRHS;
#endif
}

/* the LHS matrix for the given group */
const GlobalMatrixT& BasicSupportT::LHS(int group) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().LHS(group);
#else
#pragma unused(group)
	return *fLHS;
#endif
}

const StringT& BasicSupportT::InputFile(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return ifst;
#else
	return FEManager().InputFile();
#endif
}

ofstreamT& BasicSupportT::Output(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return *ofst;
#else
	return FEManager().Output();	
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
int BasicSupportT::RegisterOutput(const OutputSetT& output_set) const 
{
	return FEManager().RegisterOutput(output_set);
}
#endif

void BasicSupportT::WriteOutput(int ID, const dArray2DT& n_values, 
	const dArray2DT& e_values) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().WriteOutput(ID, n_values, e_values);
#else
#pragma unused(ID)
#pragma unused(n_values)
#pragma unused(e_values)
#endif
}

void BasicSupportT::WriteOutput(int ID, const dArray2DT& n_values) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().WriteOutput(ID, n_values);
#else
#pragma unused(ID)
#pragma unused(n_values)
#endif
}

/* return true if output is going to be written for the current time step */
bool BasicSupportT::WriteOutput(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return false;
#else
	return TimeManager().WriteOutput();	
#endif
}

const OutputSetT& BasicSupportT::OutputSet(int ID) const 
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#pragma unused(ID)
	ExceptionT::GeneralFail("BasicSupportT::OutputSet", "not supported");
	OutputSetT* dummy = NULL;
	return *dummy;
#else
	return FEManager().OutputSet(ID);
#endif
}

/* XDOF support */
XDOF_ManagerT& BasicSupportT::XDOF_Manager(void) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	ExceptionT::GeneralFail("BasicSupportT::XDOF_Manager", "not supported");
	XDOF_ManagerT* dummy = NULL;
	return *dummy;
#else
	return NodeManager();
#endif
}

/* the element group at the specified index in the element list */
ElementBaseT& BasicSupportT::ElementGroup(int index) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	ExceptionT::GeneralFail("BasicSupportT::ElementGroup", "not supported");
	ElementBaseT* dummy = NULL;
	return *dummy;
#else
	ElementBaseT* element = FEManager().ElementGroup(index);
	if (!element) ExceptionT::GeneralFail("BasicSupportT::ElementGroup");
	return *element;
#endif
}
