/* $Id: MFAugLagMultT.cpp,v 1.1 2004-05-06 18:55:49 cjkimme Exp $ */
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

	/* work space */
	fLHS(1, ElementMatrixT::kSymmetric),
	
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
	fNumConstraints.Dimension(numBCs);
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
				ExceptionT::BadInputValue("SetOfNodesKBCT::Initialize","Cannot get schedule %d \n",fScheduleNums[i]);
		}
		
		fNumConstraints[i] = pModel->NodeSetLength(fNodeSetIDs[i]);
		fNumConstrainedDOFs += fNumConstraints[i];
	}
	
	fNodeSets.Configure(fNumConstraints);
	
	for (int i = 0; i< numBCs; i++) {
		const iArrayT& nodeSet_i = pModel->NodeSet(fNodeSetIDs[i]);	
		fNodeSets.SetRow(i, nodeSet_i);
	}
	
	//pModel->ManyNodeSets(fNodeSetIDs, fConstraintNodes); // Get all constrained nodes
	
	int numElementGroups;
	in >> numElementGroups;
	if (numElementGroups != 1)
		ExceptionT::BadInputValue("MFAugLagMultT::EchoData","Can only handle 1 element block\n");
	
	in >> fBlockID;
	fBlockID--;
	
	in >> fk;
}

/* initialize data */
void MFAugLagMultT::Initialize(void)
{
	const char caller[] = "MFAugLagMultT::Initialize";

	/* set dimensions */
	int numDOF = rEqnos.MinorDim() + fNumConstrainedDOFs; 
	//fConstraintEqnos.Dimension(2*fNumConstrainedDOFs); // Eqnos for LaGrange Mults
	//fConstraintEqnos2D.Set(fNumConstrainedDOFs, 2, fConstraintEqnos.Pointer());
	//fFloatingDOF.Dimension(fNumConstrainedDOFs);
	//fFloatingDOF = 0;
	
	/* allocate memory for force vector */
	//fConstraintForce2D.Dimension(fNumConstrainedDOFs, 1);
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
	RaggedArray2DT<int> disp_eq, nodesInSupport;
	
	fConstraintEqnos.Dimension(2*fNumEqs);
	fConstraintEqnos2D.Set(fNumEqs, 2, fConstraintEqnos.Pointer());
	
	disp_eq.Configure(fSupportSizes, 1);
	nodesInSupport.Configure(fSupportSizes, 1);
	
	LinkedListT<int> supp_j;
	LinkedListT<double> phi_j;
	int ctr = 0;
	for (int i = 0; i < fLocallyNumberedNodeSets.MajorDim(); i++) {
		int* row_i = fLocallyNumberedNodeSets(i);
		for (int j = 0; j < fLocallyNumberedNodeSets.MinorDim(i); j++) {
		
			mfElemGroup->NodalSupportAndPhi(*row_i++, supp_j, phi_j);
			
			int* row_j = nodesInSupport(ctr++);
			supp_j.Top();
			while (supp_j.Next()) 
				*row_j++ = *(supp_j.CurrentValue());
		}
	} 
	
	fField.SetLocalEqnos(nodesInSupport, disp_eq, fConstrainedDOFs);

	/* copying eqs */
	const int* lagMultEqs = fXDOF_Nodes->XDOF_Eqnos(this, 0).Pointer();
	int* iptr = fConstraintEqnos.Pointer();
	ctr = 0;
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		int* rowptr = disp_eq(i);
		for (int j = 0; j < fNodeSets.MinorDim(i); j++) {
			for (int k = 0; k < fSupportSizes[ctr++]; k++) { 	
				*iptr++ = *rowptr++;
				*iptr++ = *lagMultEqs;
			}
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

	/* register with node manager - sets initial fConstraintDOFtags */
	iArrayT set_dims(fNumConstrainedDOFs);
	set_dims = kMFAugLagDOF;
	
	MFAugLagMultT* non_const_this = const_cast<MFAugLagMultT*>(this);
	fXDOF_Nodes->XDOF_Register(non_const_this, set_dims);	

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
	fFEManager.AssembleRHS(fGroup, fConstraintValues, fConstraintEqnos);
}

/* tangent term */
void MFAugLagMultT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused (sys_type)

	/* time integration */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer());
	
	LinkedListT<int> supp_i;
	LinkedListT<double> phi_i;
	iArrayT col_eqs(1), row_eqs(1);
 	
	double& k_entry = fLHS[0];

	/* DOF by DOF */
	const int* augLagEqnos = fXDOF_Nodes->XDOF_Eqnos(this, 0).Pointer();
	int *row_i, *row_ptr = fConstraintEqnos.Pointer();
	for (int i = 0; i < fLocallyNumberedNodeSets.MajorDim(); i++) // loop over constrained node sets
	{
		row_i = fLocallyNumberedNodeSets(i);
		
		for (int j = 0; j < fLocallyNumberedNodeSets.MinorDim(i); j++) { // loop over constrained DOFs
			
			col_eqs[0] = *augLagEqnos++;
			
			mfElemGroup->NodalSupportAndPhi(*row_i++, supp_i, phi_i);
		
			supp_i.Top(); phi_i.Top();
			while (supp_i.Next() && phi_i.Next()) {
				k_entry = *(phi_i.CurrentValue()) * constK;
				row_eqs[0] = *row_ptr++;
				row_ptr++; // have eq. nos. for LaGrange mults right now, too
			}
			
			fFEManager.AssembleLHS(fGroup, fLHS, row_eqs, col_eqs);	
			
		}
		
	}
}

void MFAugLagMultT::Reset(void)
{
	
}

void MFAugLagMultT::RegisterOutput(void) 
{
	/* initialize connectivities */
	//fContactNodes2D.Alias(fContactNodes.Length(), 1, fContactNodes.Pointer());
	
	/* output labels */
	int num_output = 1;     /* force */
	ArrayT<StringT> n_labels(num_output);
	n_labels[0] = "LAMBDA";
	
	/* register output */
	OutputSetT output_set(GeometryT::kPoint, fConstraintEqnos2D, n_labels);
	fOutputID = fFEManager.RegisterOutput(output_set);
	
}

void MFAugLagMultT::WriteOutput(ostream& out) const
{
	out << "\n M F  A u g  L a g  M u l t   D a t a :\n\n";
	
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
	LinkedListT<int> supp_j;
	LinkedListT<double> phi_j;
	int *conn_ptr = fConstraintTags.Pointer(); 
	int *tag_ptr = fConstraintDOFtags.Pointer();
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		int* row_i_global = fNodeSets(i);
		int* row_i_local = fLocallyNumberedNodeSets(i);
		for (int j = 0; j < fNodeSets.MinorDim(i); j++) {
			mfElemGroup->NodalSupportAndPhi(*row_i_local++, supp_j, phi_j);
		
			supp_j.Top();
			while (supp_j.Next()) {
				*conn_ptr++ = *(supp_j.CurrentValue());
				*conn_ptr++ = *tag_ptr;
			}
			tag_ptr++;
		}
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
 	
 	int valueIndex = 0;
 	fConstraintValues = 0.;
 	NodeManagerT* pNodeManager = fFEManager.NodeManager();
 	const ScheduleT* pSchedule;
 	for (int i = 0; i < fConstrainedDOFs.Length(); i++) {
		
		if (fCodes[i] == KBC_CardT::kFix)
			valueIndex += fNumConstraints[i];
		else {	 
			pSchedule = pNodeManager->Schedule(fScheduleNums[i]);	
			if (!pSchedule) 
				ExceptionT::GeneralFail(caller,"Cannot get schedule %d \n",fScheduleNums[i]);
			double sVal = kforce*fScales[i]*pSchedule->Value();
			for (int j = 0; j < fNumConstraints[i]; j++)
				fConstraintValues[valueIndex++] = sVal;
			
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
	
	int ctr = 0;
	fSupportSizes.Dimension(fNumConstrainedDOFs);
	for (int i = 0; i < fLocallyNumberedNodeSets.MajorDim(); i++) {
		int *row_i = fLocallyNumberedNodeSets(i);
		for (int j = 0; j < fLocallyNumberedNodeSets.MinorDim(i); j++) 
			fSupportSizes[ctr++] = mfElemGroup->SupportSize(*row_i++);
	}
	
	fNumEqs = fSupportSizes.Sum();

}