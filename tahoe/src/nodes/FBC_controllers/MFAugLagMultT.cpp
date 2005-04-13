/* $Id: MFAugLagMultT.cpp,v 1.8 2005-04-13 17:10:20 cjkimme Exp $ */
#include "MFAugLagMultT.h"

#ifdef CONTINUUM_ELEMENT

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "ModelManagerT.h"
#include "XDOF_ManagerT.h"
#include "eIntegratorT.h"
#include "FieldT.h"
#include "dMatrixT.h"
#include "ifstreamT.h"
#include "OutputSetT.h"
#include "SCNIMFT.h"

using namespace Tahoe;

/* parameters */
const int kMFAugLagDOF = 1;

/* constructor */
MFAugLagMultT::MFAugLagMultT(void):

	fNumConstrainedDOFs(0),
	
	/* RHS work space */
	fRHS_wrapper(0, fRHS),

	/* LHS work spaces */
	fLHS(ElementMatrixT::kSymmetric),
	fLHS_wrapper(0, fLHS),
	
	fOtherLHS(ElementMatrixT::kSymmetric),
	fOtherLHS_wrapper(0, fOtherLHS),
	
	fRowEqs_wrapper(0, fRowEqs),
	
	/* Pointer to the element group with access to the support of each node */
	fBlockID(-1),
	mfElemGroup(NULL)
{
	SetName("augmented_Lagrangian_KBC_meshfree");
}

/* Form of tangent matrix */
GlobalT::SystemTypeT MFAugLagMultT::TangentType(void) const
{
	return GlobalT::kSymmetric;
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

	/* the field */
	const FieldT& field = Field();

	/* dimensions */
	int ndof_u = field.NumDOF();

	/* collect displacement DOF's */
	iArrayT dofs(fNumConstrainedDOFs);
	
	fConstraintEqnos.Dimension(2*(fNumEqs+fNumConstrainedDOFs));
	fConstraintEqnos2D.Set(fNumEqs+fNumConstrainedDOFs, 2, fConstraintEqnos.Pointer());
	
	fEqNos.Configure(fSupportSizes, 1);
	
	int* constrainedDOFPtr = fConstrainedDOFs.Pointer();
	int ctr = 0;	
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		for (int j = 0; j < fNodeSets.MinorDim(i); j++)
			dofs[ctr++] = *constrainedDOFPtr;
		constrainedDOFPtr++;
	}
	
	field.SetLocalEqnos(fsupport, fEqNos, dofs);

	/* copying eqs */
	const int* lagMultEqs = FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0).Pointer();
	int* iptr = fConstraintEqnos.Pointer();
	ctr = 0;
	for (int i = 0; i < fNodeSets.MajorDim(); i++) {
		int* rowptr = fEqNos(i);
		for (int j = 0; j < fNodeSets.MinorDim(i); j++) {
			for (int k = 0; k < fSupportSizes[ctr]; k++) { 	
				*iptr++ = *rowptr++;
				*iptr++ = *lagMultEqs;
			}
		// try explicitly telling the solver about the XDOF diagonal entries
		*iptr++ = *lagMultEqs;
		*iptr++ = *lagMultEqs;
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
	constraints.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
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
	
	/* field support */
	const FieldSupportT& field_support = FieldSupport();

	/* XDOF manager */
	XDOF_ManagerT& xdof_manager = field_support.XDOF_Manager();
	
	/* assemble */
	field_support.AssembleRHS(fGroup, fConstraintValues, xdof_manager.XDOF_Eqnos(this, 0));
	
	/* get current values of constraints */
	const dArray2DT& constr = xdof_manager.XDOF(this, 0);
	dArrayT force(constr.MajorDim(), constr.Pointer());
	
	iArrayT eqNos;
	for (int i = 0; i < fNumConstrainedDOFs; i++) {
		fRHS_wrapper.SetLength(fSupportSizes[i], false);
		fRHS.Copy(fphi(i));
		fRHS *= fk*fConstraintValues[i] - constKd*force[i];
		eqNos.Set(fSupportSizes[i], fEqNos(i));
		field_support.AssembleRHS(fGroup, fRHS, eqNos);
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

	/* field support */
	const FieldSupportT& field_support = FieldSupport();

	/* XDOF manager */
	XDOF_ManagerT& xdof_manager = field_support.XDOF_Manager();

	/* DOF by DOF */
	const int* augLagEqnos = xdof_manager.XDOF_Eqnos(this, 0).Pointer();
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
		
		field_support.AssembleLHS(fGroup, fOtherLHS, fRowEqs);
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

	/* field support */
	const FieldSupportT& field_support = FieldSupport();

	/* XDOF manager */
	XDOF_ManagerT& xdof_manager = field_support.XDOF_Manager();
	
	xdof_manager.XDOF_Register(this, set_dims);	
	int max_support = fSupportSizes.Max() + 1;
	fLHS_wrapper.SetDimensions(max_support);
	fOtherLHS_wrapper.SetDimensions(max_support);
	fRowEqs_wrapper.SetLength(max_support, false);
	fRHS_wrapper.SetLength(max_support, false);
	
	/* output labels */
	int num_output = 2*Field().NumDOF();     /* force, contrained value, actual value */
	ArrayT<StringT> n_labels(num_output);
	if (num_output == 4) {
		n_labels[0] = "LAMBDA[1]";
		n_labels[1] = "LAMBDA[2]";
		n_labels[2] = "h[1]";
		n_labels[3] = "h[2]";
	} else { // 3d
		n_labels[0] = "LAMBDA[1]";
		n_labels[1] = "LAMBDA[2]";
		n_labels[2] = "LAMBDA[3]";
		n_labels[3] = "h[1]";
		n_labels[4] = "h[2]";
		n_labels[5] = "h[3]";
	}
	
	/* register output */
	fUnionOfNodes.Union(fFlattenedNodeSets);
	OutputSetT output_set(GeometryT::kPoint, fFlattenedNodeSets, n_labels);
	fOutputID = field_support.RegisterOutput(output_set);
	
	fKey.Dimension(fFlattenedNodeSets.MajorDim());
	
	// ouch! N^2 search
	for (int i = 0; i < fUnionOfNodes.Length(); i++) {
		int val_i = fUnionOfNodes[i];
		int* val_j = fFlattenedNodeSets.Pointer();
		for (int j = 0; j < fFlattenedNodeSets.MajorDim(); j++) 
			if (*val_j++ == val_i)
				fKey[j] = i;
	}
}

void MFAugLagMultT::WriteOutput(ostream& out) const
{
#pragma unused(out)
	const FieldT& field = Field();
	int ndof = field.NumDOF();
	int num_output = 2*ndof;
	// compute union of fSupportSizes 
	dArray2DT n_values(fUnionOfNodes.Length(), num_output);
	n_values = 0.;

	/* sources */
	const FieldSupportT& field_support = FieldSupport();
	XDOF_ManagerT& xdof_manager = field_support.XDOF_Manager();
	
	/* get current values of constraints */
	const dArray2DT& constr = xdof_manager.XDOF(this, 0);
	//dArray2DT us(n_values.MajorDim(), field.NumDOF());
	//mfElemGroup->InterpolatedFieldAtNodes(fLocalFlatNodes, us);
	
	int ctr = 0;
	n_values = 0.;
	double *nptr = n_values.Pointer();
	for (int i = 0; i < fNodeSetIDs.Length(); i++)	{
		int which_dof = fConstrainedDOFs[i];
		for (int j = 0; j < fLocallyNumberedNodeSets.MinorDim(i); j++, ctr++) {
			n_values(fKey[ctr],which_dof) = constr(ctr,0);
			n_values(fKey[ctr],which_dof + ndof) = fConstraintValues[ctr];
			//*nptr++ = us(ctr, which_dof);
		}
	}
	
	/* send output */
	dArray2DT e_values;
	field_support.WriteOutput(fOutputID, n_values, e_values);
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
int MFAugLagMultT::Group(void) const { return Field().Group(); };

/* describe the parameters needed by the interface */
void MFAugLagMultT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);

	list.AddParameter(fBlockID, "element_group");

	ParameterT regularization(ParameterT::Double, "regularization");
	regularization.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(regularization);
}

/* information about subordinate parameter lists */
void MFAugLagMultT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FBC_ControllerT::DefineSubs(sub_list);

	/* kinematic boundary conditions */
	sub_list.AddSub("kinematic_BC", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MFAugLagMultT::NewSub(const StringT& name) const
{
	if (name == "kinematic_BC")
	{
		/* use definition from field */
		FieldSupportT field_support;
		FieldT field(field_support);
		return field.NewSub(name);
	}
	else /* inherited */
		return FBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void MFAugLagMultT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFAugLagMultT::TakeParameterList";

	/* inherited */
	FBC_ControllerT::TakeParameterList(list);

	/* dimension */
	int numBCs = list.NumLists("kinematic_BC");
	fNodeSetIDs.Dimension(numBCs);
	fConstrainedDOFs.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	fScheduleNums.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	iArrayT fNumConstraints(numBCs);
	fScales.Dimension(numBCs);

	/* field support */
	const FieldSupportT& field_support = FieldSupport();
	ModelManagerT& model_manager = field_support.ModelManager();

	fNumConstrainedDOFs = 0;
	for (int i = 0; i < numBCs; i++) {

		/* KBC parameters */
		const ParameterListT& KBC_params = list.GetList("kinematic_BC", i);

		fNodeSetIDs[i] = KBC_params.GetParameter("node_ID");
		fConstrainedDOFs[i] = KBC_params.GetParameter("dof");
		fConstrainedDOFs[i]--;

		int typ = KBC_params.GetParameter("type");
		fCodes[i] = KBC_CardT::int2CodeT(typ + 1);
		if (fCodes[i] != KBC_CardT::kFix && fCodes[i] != KBC_CardT::kDsp)
			ExceptionT::GeneralFail(caller, "code %d must be 0 or 1", i);
	
		fScheduleNums[i] = KBC_params.GetParameter("schedule");
		fScales[i] = KBC_params.GetParameter("value");
		
		if (fCodes[i] != KBC_CardT::kFix)
		{
			fScheduleNums[i]--;
			const ScheduleT* schedule = field_support.Schedule(fScheduleNums[i]);
			if (!schedule) 
				ExceptionT::BadInputValue(caller, "cannot get schedule %d", fScheduleNums[i]);
		}
		
		fNumConstraints[i] = model_manager.NodeSetLength(fNodeSetIDs[i]);
		fNumConstrainedDOFs += fNumConstraints[i];
	}
	
	fConstraintValues.Dimension(fNumConstrainedDOFs);	
	fNodeSets.Configure(fNumConstraints);
	
	for (int i = 0; i< numBCs; i++) {
		const iArrayT& nodeSet_i = model_manager.NodeSet(fNodeSetIDs[i]);	
		fNodeSets.SetRow(i, nodeSet_i);
	}

	fBlockID = list.GetParameter("element_group");
	fBlockID--;
	fk = list.GetParameter("regularization");

	/* allocate memory for force vector */
	fConstraintForce.Dimension(fNumConstrainedDOFs);
	fConstraintForce = 0.0;
}

/**********************************************************************
 * Private                                                            *
 **********************************************************************/

void MFAugLagMultT::ComputeConstraintValues(double kforce)
{
 	const char caller[] = "MFAugLagMultT::ComputeConstraintValues";
 	
 	dArray2DT us(fLocalFlatNodes.Length(), Field().NumDOF());
	mfElemGroup->InterpolatedFieldAtNodes(fLocalFlatNodes, us);
	
 	int valueIndex = 0;
 	fConstraintValues = 0.;
 	for (int i = 0; i < fCodes.Length(); i++) {
		
		int which_dof = fConstrainedDOFs[i];
		if (fCodes[i] == KBC_CardT::kFix)
			valueIndex += fNodeSets.MinorDim(i);
		else {	 
 			const ScheduleT* schedule = FieldSupport().Schedule(fScheduleNums[i]);	
			if (!schedule) ExceptionT::GeneralFail(caller, "Cannot get schedule %d", fScheduleNums[i]);
			double sVal = fScales[i]*schedule->Value();
			for (int j = 0; j < fNodeSets.MinorDim(i); j++, valueIndex++)
				fConstraintValues[valueIndex] = -kforce*(us(valueIndex,which_dof) - sVal);
			
		}
	}
}

void MFAugLagMultT::ChatWithElementGroup(void) {

	const char caller[] = "MFAugLagMultT::ChatWithElementGroup";

	if (mfElemGroup) 
		ExceptionT::GeneralFail(caller,"MF element block communicator already exists!");

	/* communicate with the meshfree element */
	ElementBaseT& elemGroup = FieldSupport().ElementGroup(fBlockID);
	
	mfElemGroup = dynamic_cast<SCNIMFT*>(&elemGroup);
	if (!mfElemGroup)
		ExceptionT::GeneralFail(caller,"Cannot cast group %d to meshfree class",fBlockID);  
	
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

#endif /* CONTINUUM_ELEMENT */
