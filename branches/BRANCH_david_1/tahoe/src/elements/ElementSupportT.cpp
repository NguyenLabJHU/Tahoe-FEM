/* $Id: ElementSupportT.cpp,v 1.33 2004-07-15 08:25:44 paklein Exp $ */
#include "ElementSupportT.h"
#include "dArray2DT.h"

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
ElementSupportT::ElementSupportT(void)
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	ieqnos = NULL;
	iparams = NULL;
	fparams = NULL;
	fGroupAverage = new GroupAverageT;
#endif
}

/** destructor */
ElementSupportT::~ElementSupportT(void)
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	delete fGroupAverage;
#endif
}

/* return the index of the given element group */
int ElementSupportT::ElementGroupNumber(const ElementBaseT* group) const
{
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	return 0;
#else
	return FEManager().ElementGroupNumber(group);
#endif
}

#ifdef _FRACTURE_INTERFACE_LIBRARY_
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
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	nodes.ResetAverage(n_values);
#else
	fGroupAverage->ResetAverage(n_values);
#endif
}

/* assemble values */
void ElementSupportT::AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	NodeManagerT& node_man = const_cast<NodeManagerT&>(NodeManager());
	node_man.AssembleAverage(nodes, vals);
#else
    fGroupAverage->AssembleAverage(nodes, vals);
#endif
}

/* average assembled values */
const dArray2DT& ElementSupportT::OutputAverage(void) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	return nodes.OutputAverage();
#else
	return fGroupAverage->OutputAverage();
#endif
}

/* return averaged values for the nodes with assembled values */
void ElementSupportT::OutputUsedAverage(dArray2DT& average_values) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	nodes.OutputUsedAverage(average_values);
#else
	fGroupAverage->OutputUsedAverage(average_values);
#endif
}

#ifdef _FRACTURE_INTERFACE_LIBRARY_
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
#ifdef _FRACTURE_INTERFACE_LIBRARY_
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
#else
	FEManager().WriteOutput(ID, n_values, e_values);
#endif
}

/* write a snapshot */
void ElementSupportT::WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
	const dArray2DT& values, const ArrayT<StringT>& labels) const
{
	FEManager().WriteOutput(file, coords, node_map, values, labels);
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
const ArrayT<StringT>& ElementSupportT::Argv(void) const { return FEManager().Argv(); }
bool ElementSupportT::CommandLineOption(const char* str) const { return FEManager().CommandLineOption(str); }
bool ElementSupportT::CommandLineOption(const char* str, int& index) const { return FEManager().CommandLineOption(str, index); }
#endif
