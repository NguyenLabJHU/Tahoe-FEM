/* $Id: ElementSupportT.cpp,v 1.4.4.4 2002-10-18 22:37:22 cjkimme Exp $ */
#include "ElementSupportT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#ifndef _SIERRA_TEST_
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "eControllerT.h"
#include "nControllerT.h"
#include "FieldT.h"
#else
#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#endif

/* constructor */

using namespace Tahoe;

ElementSupportT::ElementSupportT(void)
{
#ifndef _SIERRA_TEST_
	/* clear */
	SetFEManager(NULL);
#else
	fNumSD = 3;
	fTimeStep = 0.;
	fItNum = 0;
	fCurrentCoordinates = NULL;
	fInitialCoordinates = NULL;
	ieqnos = NULL;
	iparams = NULL;
	fparams = new double[7];
	fparams[0] = 1.;
	fparams[1] = 0.;
	fparams[2] = fparams[3] = .05;
	fparams[4] = 1.;
	fparams[5] = 100.;
	fparams[6] = 10.;
#endif
}

#ifndef _SIERRA_TEST_
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
	else
	{	
		fNumSD = fNumNodes = 0;
	}		
}

/* Tahoe version string */
const StringT& ElementSupportT::Version(void) const
{
	return FEManager().Version();
}
#endif // _SIERRA_TEST_

bool ElementSupportT::PrintInput(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().PrintInput();
#else
	return false;
#endif
}

/*Should return something not in the NodeManager*/
const dArray2DT& ElementSupportT::InitialCoordinates(void) const
{
#ifndef _SIERRA_TEST_
	return Nodes().InitialCoordinates();
#else
	return *fInitialCoordinates;
#endif
}

const dArray2DT& ElementSupportT::CurrentCoordinates(void) const
{
#ifndef _SIERRA_TEST_
	return Nodes().CurrentCoordinates();
#else
	return *fCurrentCoordinates;
#endif
}

void ElementSupportT::RegisterCoordinates(LocalArrayT& array) const
{
#ifndef _SIERRA_TEST_
	Nodes().RegisterCoordinates(array);
#else
    switch (array.Type())
    {
    	case LocalArrayT::kInitCoords:
    	{
    		array.SetGlobal(*fInitialCoordinates);
    		break;
    	}
    	case LocalArrayT::kCurrCoords:
    	{
    		array.SetGlobal(*fCurrentCoordinates);
    		break;
    	}
    	default:
    		cout << "\n FieldSupportT::RegisterCoordinates: not a coordinate type: " 
                 << array.Type() << endl;
            throw eGeneralFail;
     }
#endif
}

/* return a  schedule function */
const ScheduleT* ElementSupportT::Schedule(int num) const
{
#ifndef _SIERRA_TEST_
	return FEManager().Schedule(num);
#else
#pragma unused(num)
	return NULL;
#endif
}

const int& ElementSupportT::IterationNumber(int group) const 
{ 
#ifndef _SIERRA_TEST_
	return FEManager().IterationNumber(group); 
#else
#pragma unused(group)
	return fItNum;
#endif
}

const char* ElementSupportT::Exception(int exception) const
{
#ifndef _SIERRA_TEST_
	return FEManager().Exception(exception);
#else
	const char* sbntma[] = {
 /* 0 */ "no error",
 /* 1 */ "general fail",
 /* 2 */ "stop",
 /* 3 */ "out of memory",
 /* 4 */ "index out of range",
 /* 5 */ "dimension mismatch",
 /* 6 */ "invalid value read from input",
 /* 7 */ "zero or negative jacobian",
 /* 8 */ "MPI message passing error",
 /* 9 */ "database read failure",
 /*10 */ "bad MP heartbeat",
 /*11 */ "unknown"}; 
	return sbntma[exception];
#endif
}

int ElementSupportT::ElementGroupNumber(const ElementBaseT* element) const
{ 
#ifndef _SIERRA_TEST_
	return FEManager().ElementGroupNumber(element); 
#else
#pragma unused(element)
	return 0;
#endif
}

const double& ElementSupportT::Time(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().Time();
#else
	return fTimeStep;
#endif
}

const double& ElementSupportT::TimeStep(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().TimeStep();
#else	
	return fTimeStep;
#endif
}

const int& ElementSupportT::StepNumber(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().StepNumber();
#else
	return fItNum;
#endif
}

const int& ElementSupportT::NumberOfSteps(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().NumberOfSteps();
#else
	return fItNum;
#endif
}

#ifndef _SIERRA_TEST_
/* the element group at the specified index in the element list */
ElementBaseT& ElementSupportT::ElementGroup(int index) const
{
	ElementBaseT* element = FEManager().ElementGroup(index);
	if (!element) throw eGeneralFail;
	return *element;
}
#endif

/* geometry information */
ModelManagerT& ElementSupportT::Model(void) const
{
#ifndef _SIERRA_TEST_
	ModelManagerT* model = FEManager().ModelManager();
	if (!model) throw eGeneralFail;
	return *model;
#else
	if (!fModelManager) throw eGeneralFail;
	return *fModelManager;
#endif
}

#ifndef _SIERRA_TEST_
/* XDOF support */
XDOF_ManagerT& ElementSupportT::XDOF_Manager(void) const
{
	return Nodes();
}
#endif

/* node number map. returns NULL if there is not a map */
const iArrayT* ElementSupportT::NodeMap(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().NodeMap();
#else
	return NULL;
#endif
}

/* return a pointer to the field */
#ifndef _SIERRA_TEST_
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

#else //_SIERRA_TEST_

void ElementSupportT::SetNumNodes(int nn)
{
	fNumNodes = nn;
}

void ElementSupportT::SetTimeStep(double dt)
{
	fTimeStep = dt;
}

/* The following two functions should be called only once to set the pointers */
void ElementSupportT::SetInitialCoordinates(dArray2DT *initialCoords)
{	
	fInitialCoordinates = initialCoords;
}

void ElementSupportT::SetCurrentCoordinates(dArray2DT* currentCoords)
{	
	fCurrentCoordinates = currentCoords;
}

/* The following two functions can be called repeatedly to change the contents of
 * the coordinate arrays.
 */
void ElementSupportT::SetInitialCoordinates(double *initialCoords)
{	
	double *finit = fInitialCoordinates->Pointer();

	for (int i = 0; i < fInitialCoordinates->Length();i++)
		*finit++ = *initialCoords++;		
}

void ElementSupportT::SetCurrentCoordinates(double *currentCoords)
{
	double *fcurr = fCurrentCoordinates->Pointer();
	
	for (int i = 0; i < fCurrentCoordinates->Length(); i++)
		*fcurr++ = *currentCoords++;
}

void ElementSupportT::UpdateCurrentCoordinates(double *displacements)
{
	double *fcurr = fCurrentCoordinates->Pointer();
	double *finit = fInitialCoordinates->Pointer();
	
	for (int i = 0; i < fCurrentCoordinates->Length(); i++)
		*fcurr++ = *finit++ + *displacements++;
}

void ElementSupportT::SetModelManager(ModelManagerT *modelManager)
{
	fModelManager = modelManager;
}

void ElementSupportT::SetNumElements(int nelem)
{
	fElem = nelem;
}	

void ElementSupportT::SetEqnos(int *conn, const int& nElem, const int& nElemNodes, 
	const int& nNodes)
{
	ieqnos = new iArrayT();
	ieqnos->Allocate(nElem*nElemNodes*3);
	int *iptr, ioff;
	iptr = ieqnos->Pointer();
	for (int i = 0; i < nElem*nElemNodes; i++)
	{
		ioff = (*conn++)*3; 
		for (int k = 0; k < 3; k++)
			*iptr++ = ioff++;
	}
	
	/* Allocate vector for residual while we're here */
	fResidual = new dArrayT();
	fResidual->Allocate(3*nNodes);
}

#endif

/* element number map for the given block ID */
const iArrayT* ElementSupportT::ElementMap(const StringT& block_ID) const
{
#ifndef _SIERRA_TEST_
	return FEManager().ElementMap(block_ID);
#else
#pragma unused(block_ID)
	return NULL;
#endif
}

#ifndef _SIERRA_TEST_
/* MP */
int ElementSupportT::Size(void) const 
{ 
#ifndef _SIERRA_TEST_
	return FEManager().Size(); 
#else
	return 1;
#endif
}

int ElementSupportT::Rank(void) const 
{
#ifndef _SIERRA_TEST_
	return FEManager().Rank();
#else
	return 1;
#endif 
}

void ElementSupportT::IncomingNodes(iArrayT& nodes_in) const
{
#ifndef _SIERRA_TEST_
	FEManager().IncomingNodes(nodes_in);
#else
#pragma unused(nodes_in)
#endif
}

void ElementSupportT::OutgoingNodes(iArrayT& nodes_out) const
{
#ifndef _SIERRA_TEST_
	FEManager().OutgoingNodes(nodes_out);
#else
#pragma unused(nodes_out)
#endif
}

void ElementSupportT::SendExternalData(const dArray2DT& all_out_data) const
{
#ifndef _SIERRA_TEST_
	FEManager().SendExternalData(all_out_data);
#else
#pragma unused(all_out_data)
#endif
}

void ElementSupportT::RecvExternalData(dArray2DT& external_data) const
{
#ifndef _SIERRA_TEST_
	FEManager().RecvExternalData(external_data);
#else
#pragma unused(external_data)
#endif
}

#endif

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& eqnos) const
{
#ifndef _SIERRA_TEST_
	FEManager().AssembleLHS(group, elMat, eqnos);
#else
#pragma unused(group)
#pragma unused(elMat)
#pragma unused(eqnos)
#endif
}

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos) const
{
#ifndef _SIERRA_TEST_
	FEManager().AssembleLHS(group, elMat, row_eqnos, col_eqnos);
#else
#pragma unused(group)
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
#endif
}

void ElementSupportT::AssembleRHS(int group, const dArrayT& elRes, 
	const nArrayT<int>& eqnos) const
{
#ifndef _SIERRA_TEST_
	FEManager().AssembleRHS(group, elRes, eqnos);
#else
#pragma unused(elRes)
#pragma unused(eqnos)
	cout <<"elRes.Length() = "<<group<<"\n";
	double *fp = elRes.Pointer();
	int *ip = ieqnos->Pointer() + group*elRes.Length();
	for (int i = 0;i < elRes.Length();i++)
		(*fResidual)[*ip++] += *fp++;
	fp = fResidual->Pointer();
	for (int i = 0;i < fResidual->Length(); i++)
		cout <<"i = "<<i<<" "<<*fp++<<"\n";
#endif
}

/* initialize work space to the number of values to be averaged */
void ElementSupportT::ResetAverage(int n_values) const
{
#ifndef _SIERRA_TEST_
	Nodes().ResetAverage(n_values);
#else
#pragma unused(n_values)
#endif
}

/* assemble values */
void ElementSupportT::AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const
{
#ifndef _SIERRA_TEST_
	Nodes().AssembleAverage(nodes, vals);
#else
#pragma unused(nodes)
#pragma unused(vals)
#endif
}

/* average assembled values */
const dArray2DT& ElementSupportT::OutputAverage(void) const
{
#ifndef _SIERRA_TEST_
	return Nodes().OutputAverage();
#else
	return *fCurrentCoordinates;
#endif
}

/* return averaged values for the nodes with assembled values */
void ElementSupportT::OutputUsedAverage(dArray2DT& average_values) const
{
#ifndef _SIERRA_TEST_
	Nodes().OutputUsedAverage(average_values);
#else
#pragma unused(average_values)
#endif
}

ifstreamT& ElementSupportT::Input(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().Input();
#else
	return *ifst;
#endif
}

ofstreamT& ElementSupportT::Output(void) const
{
#ifndef _SIERRA_TEST_
	return FEManager().Output();
#else
	return *ofst;
#endif
}

int ElementSupportT::RegisterOutput(const OutputSetT& output_set) const
{
#ifndef _SIERRA_TEST_
	return FEManager().RegisterOutput(output_set);
#else
#pragma unused(output_set)
	return 0;
#endif
}

void ElementSupportT::WriteOutput(int ID, const dArray2DT& n_values, 
	const dArray2DT& e_values) const
{
#ifndef _SIERRA_TEST_
	FEManager().WriteOutput(ID, n_values, e_values);
#else
#pragma unused(ID)
#pragma unused(n_values)
#pragma unused(e_values)
#endif
}




