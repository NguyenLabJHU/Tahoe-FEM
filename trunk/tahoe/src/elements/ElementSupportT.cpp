/* $Id: ElementSupportT.cpp,v 1.32 2004-06-26 18:29:48 paklein Exp $ */
#include "ElementSupportT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "FEManagerT.h"
#include "TimeManagerT.h"
#include "CommManagerT.h"
#include "NodeManagerT.h"
#include "eIntegratorT.h"
#include "nIntegratorT.h"
#include "FieldT.h"
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
ElementSupportT::ElementSupportT(void):
	fCurrentCoordinates(NULL),
	fInitialCoordinates(NULL)
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* clear */
	SetFEManager(NULL);
#else
	fNumSD = 3;
	fTimeStep = 0.;
	fItNum = 0;
	ieqnos = NULL;
	iparams = NULL;
	fparams = NULL;
	fGroupAverage = new GroupAverageT();
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
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

		/* set model manager */
		fModelManager = fe->ModelManager();

		/* set time manager */
		fTimeManager = fe->TimeManager();

		/* set comm manager */
		fCommManager = fe->CommManager();
	}
	else
	{
		fAnalysis = GlobalT::kNoAnalysis;
		fRunState = NULL;

		/* clear nodal information */
		SetNodes(NULL);

		/* clear model manager */
		fModelManager = NULL;
		
		/* clear time manager */
		fTimeManager = NULL;

		/* clear comm manager */
		fCommManager = NULL;
	}
}

/* (re-)set the NodeManagerT */
void ElementSupportT::SetNodes(NodeManagerT* nodes)
{
	fNodes = nodes;
	if (nodes)
	{
		fInitialCoordinates = &(nodes->InitialCoordinates());
		fCurrentCoordinates = &(nodes->CurrentCoordinates());
	}
	else
	{
		fInitialCoordinates = NULL;
		fCurrentCoordinates = NULL;
	}		
}

/* Tahoe version string */
const char* ElementSupportT::Version(void) const
{
	return FEManager().Version();
}
#endif // _FRACTURE_INTERFACE_LIBRARY_

bool ElementSupportT::PrintInput(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().PrintInput();
#else
	return false;
#endif
}

/*Should return something not in the NodeManager*/
const dArray2DT& ElementSupportT::InitialCoordinates(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return Nodes().InitialCoordinates();
#else
	return *fInitialCoordinates;
#endif
}

const dArray2DT& ElementSupportT::CurrentCoordinates(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return Nodes().CurrentCoordinates();
#else
	return *fCurrentCoordinates;
#endif
}

void ElementSupportT::RegisterCoordinates(LocalArrayT& array) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
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
            throw ExceptionT::kGeneralFail;
     }
#endif
}

/* return a  schedule function */
const ScheduleT* ElementSupportT::Schedule(int num) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Schedule(num);
#else
#pragma unused(num)
	return NULL;
#endif
}

/* return the iteration number for the current solver group */
int ElementSupportT::IterationNumber(void) const
{ 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().IterationNumber(); 
#else
	return fItNum;
#endif
}

const int& ElementSupportT::IterationNumber(int group) const 
{ 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().IterationNumber(group); 
#else
#pragma unused(group)
	return fItNum;
#endif
}

/* the group number being solved or -1 if not defined */
int ElementSupportT::CurrentGroup(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().CurrentGroup();
#else
	return -1;
#endif
}

const char* ElementSupportT::Exception(ExceptionT::CodeT exception) const
{
	return ExceptionT::ToString(exception);
}

int ElementSupportT::ElementGroupNumber(const ElementBaseT* element) const
{ 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().ElementGroupNumber(element); 
#else
#pragma unused(element)
	return 0;
#endif
}

const double& ElementSupportT::Time(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Time();
#else
	return fTimeStep;
#endif
}

const double& ElementSupportT::TimeStep(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().TimeStep();
#else	
	return fTimeStep;
#endif
}

const int& ElementSupportT::StepNumber(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().StepNumber();
#else
	return fItNum;
#endif
}

const int& ElementSupportT::NumberOfSteps(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().NumberOfSteps();
#else
	return fItNum;
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* the element group at the specified index in the element list */
ElementBaseT& ElementSupportT::ElementGroup(int index) const
{
	ElementBaseT* element = FEManager().ElementGroup(index);
	if (!element) throw ExceptionT::kGeneralFail;
	return *element;
}
#endif

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* XDOF support */
XDOF_ManagerT& ElementSupportT::XDOF_Manager(void) const
{
	return Nodes();
}
#endif

/* node number map. returns NULL if there is not a map */
const ArrayT<int>* ElementSupportT::NodeMap(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().NodeMap();
#else
	return NULL;
#endif
}

/* return a pointer to the field */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
const FieldT* ElementSupportT::Field(const char* name) const
{
	return Nodes().Field(name);
}

/* return the element controller appropriate for the given field */
const eIntegratorT* ElementSupportT::eIntegrator(const FieldT& field) const
{
	return &(field.nIntegrator().eIntegrator());
}

#else //_FRACTURE_INTERFACE_LIBRARY_

void ElementSupportT::SetNumNodes(int nn)
{
	fNumNodes = nn;
	fGroupAverage->SetNumAverageRows(fNumNodes);
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
/* Try it without copying memory. Just use set */
//    fInitialCoordinates->Set(fNumNodes,fNumSD,initialCoords);
}

void ElementSupportT::SetCurrentCoordinates(double *currentCoords)
{
	double *fcurr = fCurrentCoordinates->Pointer();
	
	for (int i = 0; i < fCurrentCoordinates->Length(); i++)
		*fcurr++ = *currentCoords++;
/* Try it without copying memory. Just use set */
//    fCurrentCoordinates->Set(fNumNodes,fNumSD,currentCoords);
}

/* This function isn't currently being used. Don't know if it needs to
 * stay around.
 */
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
#pragma unused(nNodes)
	ieqnos = new iArrayT();
	ieqnos->Dimension(nElem*nElemNodes*3);
	int *iptr, ioff;
	iptr = ieqnos->Pointer();
	for (int i = 0; i < nElem*nElemNodes; i++)
	{
		ioff = (*conn++)*fNumSD; 
		for (int k = 0; k < fNumSD; k++)
			*iptr++ = ioff++;
	}
	
	/* Allocate left- and right-hand sides while we're here */
	/* Let SIERRA control the memory for the residual */
	fResidual = new dArrayT();

#pragma message("Do I really want to allocate a stiffness matrix?")
	fStiffness = new dMatrixT(ElementMatrixT::kNonSymmetric);
//	fStiffness->Dimension(fNumSD*nNodes);
}

void ElementSupportT::SetMaterialInput(double *inputFloats, int length)
{
	fparams = new dArrayT();
	fparams->Dimension(length);
	
	double *ftmp = fparams->Pointer();
	for (int i = 0; i < length; i++)
		*ftmp++ = *inputFloats++;
	
}
	
void ElementSupportT::SetElementInput(int *inputInts, int length)
{
	iparams = new iArrayT();
	iparams->Dimension(length);
	
	int *itmp = iparams->Pointer();
	for (int i = 0; i < length; i++)
		*itmp++ = *inputInts++;
}

int ElementSupportT::ReturnInputInt(CodeT label) 
{ 
		return (*iparams)[label];
}

void ElementSupportT::SetResidual(double *nodalForces)
{
	fResidual->Set(fNumSD*fNumNodes,nodalForces);
}

void ElementSupportT::SetStateVariableArray(double *incomingArray)
{
	fStateVars = incomingArray;
}

double *ElementSupportT::StateVariableArray(void)
{
	return fStateVars;
}

void ElementSupportT::SetBlockID(StringT& Id)
{
	sBlockID = Id;
}

StringT& ElementSupportT::BlockID(void)
{
	return sBlockID;
}

void ElementSupportT::OutputSize(int& nNodeOutputVars, int& nElemOutputVars)
{
	nNodeOutputVars = fNodeOutputLabels.Length();
	nElemOutputVars = fElemOutputLabels.Length();
}
	
void ElementSupportT::SetOutputCodes(iArrayT& fNodalOutputCodes, iArrayT& fElementOutputCodes)
{
#pragma message("Must read in IO codes somehow")
	fNodalOutputCodes = IOBaseT::kAtInc;
	fElementOutputCodes = IOBaseT::kAtInc;
}

void ElementSupportT::SetOutputPointers(double *nodalOutput, double *elemOutput)
{
	fNodalOutput = nodalOutput;
	fElemOutput = elemOutput;
}

#endif

/* element number map for the given block ID */
const iArrayT* ElementSupportT::ElementMap(const StringT& block_ID) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().ElementMap(block_ID);
#else
#pragma unused(block_ID)
	return NULL;
#endif
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* MP */
int ElementSupportT::Size(void) const 
{ 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Size(); 
#else
	return 1;
#endif
}

int ElementSupportT::Rank(void) const 
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Rank();
#else
	return 0;
#endif 
}

/* low-level communicator */
const CommunicatorT& ElementSupportT::Communicator(void) const
{
	if (!fCommManager) 
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		ExceptionT::GeneralFail("ElementSupportT::Communicator", "pointer not set");
#else
		ExceptionT::GeneralFail("ElementSupportT::Communicator", "not supported");
#endif

	return fCommManager->Communicator();
}

const ArrayT<int>* ElementSupportT::ExternalNodes(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	if (fCommManager)
		return fCommManager->ExternalNodes();
	else
		return NULL;
#else
	return NULL;
#endif
}

const ArrayT<int>* ElementSupportT::BorderNodes(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	if (fCommManager)
		return fCommManager->BorderNodes();
	else
		return NULL;
#else
	return NULL;
#endif
}
#endif

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
	const nArrayT<int>& eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleLHS(group, elMat, eqnos);
#else
#pragma unused(eqnos)
#pragma message("ElementSupportT::AssembleLHS only fullMatrix so far")
/* NB that group is really the element number; it's an offset in my eq array */
	double *fp = elMat.Pointer();
	int *ip1 = ieqnos->Pointer() + group*elMat.Rows();
	int nElemDOF = elMat.Rows();
	
	for (int i = 0;i < elMat.Rows();i++)
	{
		int *ip2 = ieqnos->Pointer() + nElemDOF*group;

		/* go to right row of stiffness matrix */
		double *fstiffptr = fStiffness->Pointer()+(*ip1++)*fStiffness->Rows();
		*(fstiffptr + *ip2++) += *fp++;
	}
	fp = fStiffness->Pointer();
	for (int i = 0;i < fStiffness->Length(); i++)
		cout <<"i = "<<i<<" "<<*fp++<<"\n";
#endif
}

void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, 
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

void ElementSupportT::AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, 
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

void ElementSupportT::AssembleRHS(int group, const nArrayT<double>& elRes, 
	const nArrayT<int>& eqnos) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().AssembleRHS(group, elRes, eqnos);
#else
#pragma unused(eqnos)
/* NB that group is really the element number; it's an offset in my eq array */
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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	Nodes().ResetAverage(n_values);
#else
	fGroupAverage->ResetAverage(n_values);
#endif
}

/* assemble values */
void ElementSupportT::AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	Nodes().AssembleAverage(nodes, vals);
#else
    fGroupAverage->AssembleAverage(nodes,vals);
#endif
}

/* average assembled values */
const dArray2DT& ElementSupportT::OutputAverage(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return Nodes().OutputAverage();
#else
	return fGroupAverage->OutputAverage();
#endif
}

/* return averaged values for the nodes with assembled values */
void ElementSupportT::OutputUsedAverage(dArray2DT& average_values) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	Nodes().OutputUsedAverage(average_values);
#else
	fGroupAverage->OutputUsedAverage(average_values);
#endif
}

ifstreamT& ElementSupportT::Input(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Input();
#else
	return *ifst;
#endif
}

ofstreamT& ElementSupportT::Output(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return FEManager().Output();
#else
	return *ofst;
#endif
}

/* format of the output files */
IOBaseT::FileTypeT ElementSupportT::OutputFormat(void) const { return FEManager().OutputFormat(); }

#ifndef _FRACTURE_INTERFACE_LIBRARY_
int ElementSupportT::RegisterOutput(const OutputSetT& output_set) const
{
	return FEManager().RegisterOutput(output_set);
}
#else
int ElementSupportT::RegisterOutput(ArrayT<StringT>& n_labels, 
	ArrayT<StringT>& e_labels)
{
	/* copy labels */
	fNodeOutputLabels.Dimension(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = n_labels[i];
	fElemOutputLabels.Dimension(e_labels.Length());
	for (int i = 0; i < fElemOutputLabels.Length(); i++)
		fElemOutputLabels[i] = e_labels[i];
		
	return 0;
}
#endif

void ElementSupportT::WriteOutput(int ID, const dArray2DT& n_values, 
	const dArray2DT& e_values) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	FEManager().WriteOutput(ID, n_values, e_values);
#else
#pragma unused(ID)
	double *ftmp1, *ftmp2;
	ftmp1 = fNodalOutput;
	ftmp2 = n_values.Pointer();
	for (int i = 0; i < n_values.Length(); i++)
		*ftmp1++ = *ftmp2++;
	ftmp1 = fElemOutput;
	ftmp2 = e_values.Pointer();
	for (int i = 0; i < e_values.Length(); i++)
		*ftmp1++ = *ftmp2++;
#endif
}

/* return true if output is going to be written for the current time step */
bool ElementSupportT::WriteOutput(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	return TimeManager().WriteOutput();
#else
	return false;
#endif
}

/* write a snapshot */
void ElementSupportT::WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
	const dArray2DT& values, const ArrayT<StringT>& labels) const
{
	FEManager().WriteOutput(file, coords, node_map, values, labels);
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
const OutputSetT& ElementSupportT::OutputSet(int ID) const
{
	return FEManager().OutputSet(ID);
}
#endif

#ifndef _FRACTURE_INTERFACE_LIBRARY_
const ArrayT<StringT>& ElementSupportT::Argv(void) const { return FEManager().Argv(); }
bool ElementSupportT::CommandLineOption(const char* str) const { return FEManager().CommandLineOption(str); }
bool ElementSupportT::CommandLineOption(const char* str, int& index) const { return FEManager().CommandLineOption(str, index); }
#endif