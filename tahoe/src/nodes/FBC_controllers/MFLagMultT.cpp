/* $Id: MFLagMultT.cpp,v 1.1 2004-05-06 18:54:47 cjkimme Exp $ */
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

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
MFAugLagMultT::MFAugLagMultT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes,
	const FieldT& field, const dArray2DT& coords, const dArray2DT& disp):
	FBC_ControllerT(fe_manager, 
		field.Group()),

	/* references to NodeManagerT data */
	rEqnos(field.Equations()),
	rCoords(coords),
	rDisp(disp),
	fXDOF_Nodes(XDOF_nodes),
	fField(field),

	fNumConstrainedDOFs(0),

	/* work space */
	fLHS(1, ElementMatrixT::kSymmetric),
	
	mfElemGroup(NULL)
{
	
}

/* get input data */
void MFAugLagMultT::EchoData(ifstreamT& in, ostream& out)
{
#pragma unused(out)

	int numBCs;
	in >> numBCs;
	if (numBCs <= 0)
		ExceptionT::BadInputValue("MFAugLagMultT::Initialize","Expecting numBCs > 0");
	
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
		
		fNumConstraints[i] = pModel->NodeSet(fNodeSetIDs[i]).Length();
		fNumConstrainedDOFs += fNumConstraints[i];
	}
	
	// convert from nodes to DOFs
	
	int numElementGroups;
	in >> numElementGroups;
	if (numElementGroups != 1)
		ExceptionT::BadInputValue("MFAugLagMultT::EchoData","Can only handle 1 element block\n");
	
	in >> fBlockID;
}

/* initialize data */
void MFAugLagMultT::Initialize(void)
{
	const char caller[] = "MFAugLagMultT::Initialize";

	/* set dimensions */
	int numDOF = rEqnos.MinorDim() + fNumConstrainedDOFs; 
	fConstraintEqnos.Dimension(fNumConstrainedDOFs); // Eqnos for LaGrange Mults
	fConstraintEqnos2D.Set(fNumConstrainedDOFs, 1, fConstraintEqnos.Pointer());
	fFloatingDOF.Dimension(fNumConstrainedDOFs);
	fFloatingDOF = 0;
	
	/* allocate memory for force vector */
	fConstraintForce2D.Dimension(fNumConstrainedDOFs, 1);
	fConstraintForce.Set(fNumConstrainedDOFs, fConstraintForce2D.Pointer());
	fConstraintForce2D = 0.0;

	/* register with node manager - sets initial fConstraintDOFtags */
	iArrayT set_dims(fNumConstrainedDOFs);
	set_dims = kNumAugLagDOF;
	fXDOF_Nodes->XDOF_Register(this, set_dims);	
	
	/* communicate with the meshfree element */
	ElementBaseT* elemGroup = fFEManager.ElementGroup(fBlockID);
	if (!elemGroup) 
		ExceptionT::GeneralFail(caller,"Element group %d does not exist\n", fBlockID);
	
	mfElemGroup = dynamic_cast<SCNIMFT*>(elemGroup);
	if (!mfElemGroup)
		ExceptionT::GeneralFail(caller,"Cannot cast group %d to meshfree class\n",fBlockID);  
	
	mfElemGroup->GlobalToLocalNumbering(fConstraintNodes);

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
#pragma message("Rewrite MFAugLagMultT::Equations")
	/* dimensions */
	int ndof_u = rCoords.MinorDim();

	/* collect displacement DOF's */
	iArray2DT disp_eq(fConstraintNodes.Length(), ndof_u);
	fField.SetLocalEqnos(fConstraintNodes, disp_eq);

	int eq_col = 0;
	iArrayT eq_temp(fConstraintNodes.Length());

	/* displacement equations */
	fFloatingDOF = 0;
	for (int i = 0; i < ndof_u; i++)
	{
		disp_eq.ColumnCopy(i, eq_temp);
		fConstraintEqnos2D.SetColumn(eq_col++, eq_temp);

		/* check for floating DOF's */
		for (int j = 0; j < eq_temp.Length(); j++)
			if (eq_temp[j] < 1)
				fFloatingDOF[j] = 1; /* mark */
	}

	/* warning */
	if (fFloatingDOF.HasValue(1))
		cout << "\n MFAugLagMultT::Equations: node with constraint has prescribed DOF\n" 
		     <<   "     Stiffness may be approximate." << endl;	

	/* constraint equations */
	const iArray2DT& auglageqs = fXDOF_Nodes->XDOF_Eqnos(this, 0);
	for (int j = 0; j < auglageqs.MinorDim(); j++)
	{
		auglageqs.ColumnCopy(j, eq_temp);
		fConstraintEqnos2D.SetColumn(eq_col++, eq_temp);
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

void MFAugLagMultT::CloseStep(void)
{
	/* inherited */
	FBC_ControllerT::CloseStep();

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
	
	LinkedListT<int>* supp_i;
	LinkedListT<double>* phi_i;
	iArrayT col_eqs(1), row_eqs(1);
 	
	double& k_entry = fLHS[0];

	/* DOF by DOF */
	for (int i = 0; i < fNumConstrainedDOFs; i++)
	{
		col_eqs[0] = 0; 
		
		supp_i = NULL;
		phi_i = NULL;
		mfElemGroup->NodalSupportAndPhi(fConstraintNodes[i], supp_i, phi_i);
		
		if ((supp_i == NULL) || (phi_i == NULL))
			ExceptionT::GeneralFail("MFAugLagMultT::ApplyLHS","Cannot get shape function information from element block\n");
		
		supp_i->Top(); phi_i->Top();
		while (supp_i->Next() && phi_i->Next()) {
			k_entry = *(phi_i->CurrentValue()) * constK;
			row_eqs[0] = 0;
		}
			
		fFEManager.AssembleLHS(fGroup, fLHS, row_eqs, col_eqs);	
	}
}

/* returns the array for the DOF tags needed for the current config */
void MFAugLagMultT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fConstraintDOFtags.Dimension(fConstrainedDOFs.Length());
}

iArrayT& MFAugLagMultT::DOFTags(int tag_set)
{
#pragma unused (tag_set)
	return fConstraintDOFtags;
}

/* generate nodal connectivities - does nothing here */
void MFAugLagMultT::GenerateElementData(void)
{
	/* allocate space */
	fConstraintTags.Dimension(fConstrainedDOFs.Length(), 2);
	
	/* collect tags - {contact node, DOF tag} */
	fConstraintTags.SetColumn(0, fConstrainedDOFs);
	fConstraintTags.SetColumn(1, fConstraintDOFtags);
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