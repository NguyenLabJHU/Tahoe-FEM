/* $Id: MFAugLagMultT.cpp,v 1.3 2004-05-14 23:07:34 cjkimme Exp $ */
#include "MFAugLagMultT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "ModelManagerT.h"
#include "XDOF_ManagerT.h"
#include "eIntegratorT.h"
#include "FieldT.h"
#include "SCNIMFT.h"
#include "dMatrixT.h"
#include "ifstreamT.h"
#include "OutputSetT.h"

using namespace Tahoe;

/* parameters */
const int kMFAugLagDOF = 1;

/* constructor */
MFAugLagMultT::MFAugLagMultT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes,
	const FieldT& field, const dArray2DT& coords, const dArray2DT& disp):
	FBC_ControllerT(fe_manager, 
		field.Group()),

	/* references to NodeManagerT data */
	rEqnos(field.Equations()),
	rDisp(disp),
	fXDOF_Nodes(XDOF_nodes),
	fField(field),

	fNumConstrainedDOFs(0),
	
	/* RHS work space */
	fRHS(),
	fRHS_wrapper(0, fRHS),

	/* LHS work spaces */
	fLHS(ElementMatrixT::kSymmetric),
	fLHS_wrapper(0, fLHS),
	
	fOtherLHS(ElementMatrixT::kSymmetric),
	fOtherLHS_wrapper(0, fOtherLHS),
	
	fRowEqs(),
	fRowEqs_wrapper(0, fRowEqs),
	
	/* Pointer to the element group with access to the support of each node */
	mfElemGroup(NULL)
{
#pragma unused(coords)
	
}

/* Form of tangent matrix */
GlobalT::SystemTypeT MFAugLagMultT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* get input data */
void MFAugLagMultT::EchoData(ifstreamT& in, ostream& out)
{
#pragma unused(out)

	int numBCs;
	in >> numBCs;
	if (numBCs <= 0)
		ExceptionT::BadInputValue("MFAugLagMultT::Initialize","Expecting numBCs > 0");
	
	fNodeSetIDs.Dimension(numBCs);
	fConstrainedDOFs.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	fScheduleNums.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	iArrayT fNumConstraints(numBCs);
	fScales.Dimension(numBCs);
	
	const ScheduleT* pSchedule;
	NodeManagerT* pNodeManager = fFEManager.NodeManager();
	ModelManagerT* pModel = fFEManager.ModelManager();  
	  
	fNumConstrainedDOFs = 0;
	for (int i = 0; i < numBCs; i++) {
		in >> fNodeSetIDs[i];
	
		in >> fConstrainedDOFs[i];
		fConstrainedDOFs[i]--;
	
		in >> fCodes[i];
		
		if (fCodes[i] != 0 && fCodes[i] != 1)
			ExceptionT::GeneralFail("MFAugLagMultT::EchoData","Code %d must be 0 or 1\n",i);
	
		in >> fScheduleNums[i] >> fScales[i]; 
		
		if (fCodes[i] != KBC_CardT::kFix)
		{
			fScheduleNums[i]--;
			pSchedule = pNodeManager->Schedule(fScheduleNums[i]);	
			if (!pSchedule) 
				ExceptionT::BadInputValue("MFAugLagMultT::EchoData","Cannot get schedule %d \n",fScheduleNums[i]);
		}
		
		fNumConstraints[i] = pModel->NodeSetLength(fNodeSetIDs[i]);
		fNumConstrainedDOFs += fNumConstraints[i];
	}
	
	fConstraintValues.Dimension(fNumConstrainedDOFs);
	
	fNodeSets.Configure(fNumConstraints);
	
	for (int i = 0; i< numBCs; i++) {
		const iArrayT& nodeSet_i = pModel->NodeSet(fNodeSetIDs[i]);	
		fNodeSets.SetRow(i, nodeSet_i);
	}
	
	int numElementGroups;
	in >> numElementGroups;
	if (numElementGroups != 1)
		ExceptionT::BadInputValue("MFAugLagMultT::EchoData","Can only handle 1 element block\n");
	
	in >> fBlockID;
	fBlockID--;
	
	in >> fk;
	if (fk < 0.) ExceptionT::GeneralFail("MFAugLagMultT::EchoData","fk must be non-negative\n");
}

/* initialize data */
void MFAugLagMultT::Initialize(void)
{
	const char caller[] = "MFAugLagMultT::Initialize";
	
	/* allocate memory for force vector */
	fConstraintForce.Dimension(fNumConstrainedDOFs);
	fConstraintForce = 0.0;

}

void MFAugLagMultT::SetEquationNumbers(void)
{
// don't need to set FBC destinations as with the base class because
// the class collects the global equation numbers during SendEqnsToSolver()
}

/* append element equations numbers to the list */
void MFAugLagMultT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* dimensions */
	int ndof_u = rDisp.MinorDim();

	/* collect displacement DOF's */
	iArrayT dofs(fNumConstrainedDOFs);
	
	fConstraintEqnos.Dimension(2*fNumEqs);
	fConstraintEqnos2D.Set(fNumEqs, 2, fConstraintEqnos.Pointer());
	
	fEqNos.Configure(fSupportSizes, 1);
	
	int* constrainedDOFPtr = fConstrainedDOFs.Pointer();
	int ctr = 0;	
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		for (int j = 0; j < fNodeSets.MinorDim(i); j++)
			dofs[ctr++] = *constrainedDOFPtr;
		constrainedDOFPtr++;
	}
	
	fField.SetLocalEqnos(fsupport, fEqNos, dofs);

	/* copying eqs */
	const int* lagMultEqs = fXDOF_Nodes->XDOF_Eqnos(this, 0).Pointer();
	int* iptr = fConstraintEqnos.Pointer();
	ctr = 0;
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		int* rowptr = fEqNos(i);
		for (int j = 0; j < fNodeSets.MinorDim(i); j++) {
			for (int k = 0; k < fSupportSizes[ctr]; k++) { 	
				*iptr++ = *rowptr++;
				*iptr++ = *lagMultEqs;
			}
		ctr++;
		lagMultEqs++;
		}
	}
	
	/* send to solver */
	eq_1.Append(&fConstraintEqnos2D);
}

void MFAugLagMultT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)

	connects_1.Append(&fConstraintTags);
}

void MFAugLagMultT::ReadRestart(istream& in)
{
	/* inherited */
	FBC_ControllerT::ReadRestart(in);

	//in >> fLastDOF; // previous solution
}

void MFAugLagMultT::WriteRestart(ostream& out) const
{
	/* inherited */
	FBC_ControllerT::WriteRestart(out);

	//out << fLastDOF; // previous solution
}

void MFAugLagMultT::InitStep(void)
{
}

void MFAugLagMultT::CloseStep(void)
{

	/* store last converged DOF array */
	dArrayT constraints;
	constraints.Alias(fXDOF_Nodes->XDOF(this, 0));
	fLastDOF = constraints;
}

/* restore the DOF values to the last converged solution */
void MFAugLagMultT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused (tag_set)

	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

void MFAugLagMultT::InitialCondition(void)
{
	
}

/* compute the nodal contribution to the residual force vector */
void MFAugLagMultT::ApplyRHS(void)
{
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* recompute contstraint forces */
	ComputeConstraintValues(constKd);
	
	/* assemble */
	fFEManager.AssembleRHS(fGroup, fConstraintValues, fXDOF_Nodes->XDOF_Eqnos(this, 0));
	
	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this, 0);
	dArrayT force(constr.MajorDim(), constr.Pointer());
	
	iArrayT eqNos;
	for (int i = 0; i < fNumConstrainedDOFs; i++) {
		fRHS_wrapper.SetLength(fSupportSizes[i], false);
		fRHS.Copy(fphi(i));
		fRHS *= fk*fConstraintValues[i] - constKd*force[i];
		eqNos.Set(fSupportSizes[i], fEqNos(i));
		fFEManager.AssembleRHS(fGroup, fRHS, eqNos);
	}
	
}

/* tangent term */
void MFAugLagMultT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused (sys_type)

	/* time integration */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* DOF by DOF */
	const int* augLagEqnos = fXDOF_Nodes->XDOF_Eqnos(this, 0).Pointer();
	for (int i = 0; i < fNumConstrainedDOFs; i++) 
	{
		int supp = fSupportSizes[i];
					
		fLHS_wrapper.SetDimensions(supp);
		fOtherLHS_wrapper.SetDimensions(supp+1);
		fRowEqs_wrapper.SetLength(supp, false);
		
		fOtherLHS = 0.;
		fLHS.Outer(fphi(i),fphi(i),fk);
		fOtherLHS.SetBlock(0,0,fLHS);
		fLHS_wrapper.SetDimensions(supp,1);
		fLHS.Copy(fphi(i));
		fOtherLHS.SetBlock(0,supp,fLHS);
		fLHS_wrapper.SetDimensions(1,supp);
		fOtherLHS.SetBlock(supp,0,fLHS);
		fOtherLHS *= constK;
		fRowEqs.Copy(fEqNos(i));
		fRowEqs_wrapper.SetLength(supp+1, false);
		fRowEqs.Last() = augLagEqnos[i];
		
		fFEManager.AssembleLHS(fGroup, fOtherLHS, fRowEqs);
	}
}

void MFAugLagMultT::Reset(void)
{
	
}

void MFAugLagMultT::RegisterOutput(void) 
{
	/* initialize connectivities */
	//fContactNodes2D.Alias(fContactNodes.Length(), 1, fContactNodes.Pointer());
	/* register with node manager - sets initial fConstraintDOFtags */
	iArrayT set_dims(1);
	set_dims = kMFAugLagDOF;
	
	fXDOF_Nodes->XDOF_Register(this, set_dims);	
	
	int max_support = fSupportSizes.Max() + 1;
	fLHS_wrapper.SetDimensions(max_support);
	fOtherLHS_wrapper.SetDimensions(max_support);
	fRowEqs_wrapper.SetLength(max_support, false);
	fRHS_wrapper.SetLength(max_support, false);
	
	/* output labels */
	int num_output = 2;     /* force, contrained value, actual value */
	ArrayT<StringT> n_labels(num_output);
	n_labels[0] = "LAMBDA";
	n_labels[1] = "h";
	
	/* register output */
	OutputSetT output_set(GeometryT::kPoint, fFlattenedNodeSets, n_labels);
	fOutputID = fFEManager.RegisterOutput(output_set);
	
}

void MFAugLagMultT::WriteOutput(ostream& out) const
{
	out << "\n M F  A u g  L a g  M u l t   D a t a :\n\n";
	
	int num_output = 2;
	dArray2DT n_values(fSupportSizes.Length(), num_output);
	n_values = 0.;
	
	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this, 0);
	dArray2DT us(n_values.MajorDim(), rDisp.MinorDim());
	mfElemGroup->InterpolatedFieldAtNodes(fLocalFlatNodes, us);
	
	int ctr = 0;
	double *nptr = n_values.Pointer();
	for (int i = 0; i < fNodeSetIDs.Length(); i++)	{
		int which_dof = fConstrainedDOFs[i];
		for (int j = 0; j < fLocallyNumberedNodeSets.MinorDim(i); j++, ctr++) {
			*nptr++ = constr(ctr,0);
			*nptr++ = fConstraintValues[ctr];
			//*nptr++ = us(ctr, which_dof);
		}
	}
	
	/* send output */
	dArray2DT e_values;
	fFEManager.WriteOutput(fOutputID, n_values, e_values);
}

/* returns the array for the DOF tags needed for the current config */
void MFAugLagMultT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fConstraintDOFtags.Dimension(fNumConstrainedDOFs);
}

iArrayT& MFAugLagMultT::DOFTags(int tag_set)
{
#pragma unused (tag_set)
	return fConstraintDOFtags;
}

/* generate nodal connectivities  */
void MFAugLagMultT::GenerateElementData(void)
{
	ChatWithElementGroup();

	/* allocate space */
	fConstraintTags.Dimension(fNumEqs, 2);
	
	/* collect tags - {node in support of constrained node, DOF tag} */
	int *conn_ptr = fConstraintTags.Pointer(); 
	int *tag_ptr = fConstraintDOFtags.Pointer();
	for (int i = 0; i < fsupport.MajorDim(); i++) {
		int *supp_i = fsupport(i);
		for (int j = 0; j < fsupport.MinorDim(i); j++) {
				*conn_ptr++ = *supp_i++;
				*conn_ptr++ = *tag_ptr;
		}
		tag_ptr++;
	}
}

/* return the contact elements */
const iArray2DT& MFAugLagMultT::DOFConnects(int tag_set) const
{
#pragma unused (tag_set)
	return fConstraintTags;
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int MFAugLagMultT::Reconfigure(void) { return 0; }

/* return the equation group */
int MFAugLagMultT::Group(void) const { return fField.Group(); };

/**********************************************************************
 * Private                                                            *
 **********************************************************************/

void MFAugLagMultT::ComputeConstraintValues(double kforce)
{
 	const char caller[] = "MFAugLagMultT::ComputeConstraintValues\n";
 	
 	dArray2DT us(fLocalFlatNodes.Length(), rDisp.MinorDim());
	mfElemGroup->InterpolatedFieldAtNodes(fLocalFlatNodes, us);
	
 	int valueIndex = 0;
 	fConstraintValues = 0.;
 	NodeManagerT* pNodeManager = fFEManager.NodeManager();
 	const ScheduleT* pSchedule;
 	for (int i = 0; i < fCodes.Length(); i++) {
		
		int which_dof = fConstrainedDOFs[i];
		if (fCodes[i] == KBC_CardT::kFix)
			valueIndex += fNodeSets.MinorDim(i);
		else {	 
			pSchedule = pNodeManager->Schedule(fScheduleNums[i]);	
			if (!pSchedule) 
				ExceptionT::GeneralFail(caller,"Cannot get schedule %d \n",fScheduleNums[i]);
			double sVal = fScales[i]*pSchedule->Value();
			for (int j = 0; j < fNodeSets.MinorDim(i); j++, valueIndex++)
				fConstraintValues[valueIndex] = -kforce*(us(valueIndex,which_dof) - sVal);
			
		}
	}
	
}

void MFAugLagMultT::ChatWithElementGroup(void) {

	const char caller[] = "MFAugLagMultT::ChatWithElementGroup";

	if (mfElemGroup) 
		ExceptionT::GeneralFail(caller,"MF element block communicator already exists!\n");

	/* communicate with the meshfree element */
	ElementBaseT* elemGroup = fFEManager.ElementGroup(fBlockID);
	if (!elemGroup) 
		ExceptionT::GeneralFail(caller,"Element group %d does not exist\n", fBlockID);
	
	mfElemGroup = dynamic_cast<SCNIMFT*>(elemGroup);
	if (!mfElemGroup)
		ExceptionT::GeneralFail(caller,"Cannot cast group %d to meshfree class\n",fBlockID);  
	
	fLocallyNumberedNodeSets = fNodeSets; // ouch! two of these data structures
	mfElemGroup->GlobalToLocalNumbering(fLocallyNumberedNodeSets);
	
	fFlattenedNodeSets.Dimension(fNumConstrainedDOFs, 1);
	fLocalFlatNodes.Dimension(fNumConstrainedDOFs);
	int ctr = 0;
	for (int i = 0; i < fNodeSets.MajorDim(); i++)
		for (int j = 0; j < fNodeSets.MinorDim(i); j++) {
			fFlattenedNodeSets[ctr] = fNodeSets(i,j);
			fLocalFlatNodes[ctr++] = fLocallyNumberedNodeSets(i,j);
		}
		
	mfElemGroup->NodalSupportAndPhi(fLocalFlatNodes, fsupport, fphi);
	
	fSupportSizes.Dimension(fNumConstrainedDOFs);
	for (int i = 0; i < fNumConstrainedDOFs; i++)
		fSupportSizes[i] = fsupport.MinorDim(i);
			
	fNumEqs = fSupportSizes.Sum();

}