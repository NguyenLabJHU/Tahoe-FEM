/* $Id: FieldMFAugLagMultT.cpp,v 1.2 2005-06-30 21:48:23 jcmach Exp $ */
#include "FieldMFAugLagMultT.h"

#ifdef CONTINUUM_ELEMENT

#include "ParameterUtils.h"
#include "FieldSupportT.h"
#include "ModelManagerT.h"
#include "FieldT.h"
#include "SCNIMFT.h"
#include "ScheduleT.h"

using namespace Tahoe;

/* constructor */
FieldMFAugLagMultT::FieldMFAugLagMultT(void):
	fSchedule(NULL)
{
	SetName("field_augmented_Lagrangian_KBC_meshfree");
}

/* describe the parameters needed by the interface */
void FieldMFAugLagMultT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	MFAugLagMultT::DefineParameters(list);

	/* solution schedule */
	list.AddParameter(ParameterT::Integer, "solution_schedule");
}

/* information about subordinate parameter lists */
void FieldMFAugLagMultT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MFAugLagMultT::DefineSubs(sub_list);

	/* remove inherited kinematic boundary conditions */
	if (!sub_list.RemoveSub("kinematic_BC"))
		ExceptionT::GeneralFail("FieldMFAugLagMultT::DefineSubs",
			"could not remove sub \"kinematic_BC::DefineSubs\"");

	/* list of affected nodes */
	sub_list.AddSub("prescribed_node_ID_list");
}

/* accept parameter list */
void FieldMFAugLagMultT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FieldMFAugLagMultT::TakeParameterList";

	/* inherited */
	MFAugLagMultT::TakeParameterList(list);

	/* field support */
	const FieldSupportT& field_support = FieldSupport();
	ModelManagerT& model_manager = field_support.ModelManager();

	/* get nodes */
	ArrayT<StringT> id_list;
	StringListT::Extract(list.GetList("prescribed_node_ID_list"), id_list);
	model_manager.ManyNodeSets(id_list, fAllNodes);

	/* register the collective node set */
	StringT node_ID;
	if (id_list.Length() > 1) {
		node_ID = model_manager.FreeNodeSetID("80");
		if (!model_manager.RegisterNodeSet(node_ID, fAllNodes, false))
			ExceptionT::GeneralFail(caller, "could not register node set \"%s\"",
				node_ID.Pointer());
	}
	else
		node_ID = id_list[0];

	/* dimension */
	int ndof = Field().NumDOF();
	int numBCs = (id_list.Length() > 0) ? ndof : 0; /* base class assumed one KCB_CardT per BC */
	fNodeSetIDs.Dimension(numBCs);
	fConstrainedDOFs.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	fScheduleNums.Dimension(numBCs);
	fCodes.Dimension(numBCs);
	iArrayT num_constraints(numBCs);
	fScales.Dimension(numBCs);
	fPrescribedField.Dimension(fAllNodes.Length(), ndof);

	/* schedule to scale the solution */
	int schedule = list.GetParameter("solution_schedule");
	schedule--;
	fSchedule = FieldSupport().Schedule(schedule);
	if (!fSchedule) ExceptionT::GeneralFail(caller, "could not resolve schedule %d", schedule);

	fNumConstrainedDOFs = 0;
	int id_num = 0;
	for (int i = 0; i < numBCs; i++) {

		/* KBC parameters */
		fNodeSetIDs[i] = node_ID;
		fConstrainedDOFs[i] = i;
		fCodes[i] = KBC_CardT::kDsp;	
		fScheduleNums[i] = schedule;
		fScales[i] = 0.0; /* not used */

		/* count constraints */
		num_constraints[i] = model_manager.NodeSetLength(fNodeSetIDs[i]);
		fNumConstrainedDOFs += num_constraints[i];
	}
	
	fConstraintValues.Dimension(fNumConstrainedDOFs);	
	fNodeSets.Configure(num_constraints);
	
	for (int i = 0; i< numBCs; i++) {
		const iArrayT& nodeSet_i = model_manager.NodeSet(fNodeSetIDs[i]);	
		fNodeSets.SetRow(i, nodeSet_i);
	}

	/* allocate memory for force vector */
	fConstraintForce.Dimension(fNumConstrainedDOFs);
	fConstraintForce = 0.0;	
}

/**********************************************************************
 * Private                                                            *
 **********************************************************************/

void FieldMFAugLagMultT::ComputeConstraintValues(double kforce)
{
 	const char caller[] = "FieldMFAugLagMultT::ComputeConstraintValues";
 	
	/* compute field at nodes */
	int ndof = Field().NumDOF();
 	dArray2DT us(fAllNodes.Length(), ndof);
	mfElemGroup->InterpolatedFieldAtNodes(fAllNodes, us);

	/* compute the prescribed field values */
	const dArray2DT& initial_coordinates = FieldSupport().InitialCoordinates();
	for (int i = 0; i < fAllNodes.Length(); i++) {
		double scale = fSchedule->Value();
		int nd = fAllNodes[i];
// 		fPrescribedField(i,0) = scale*0.010*initial_coordinates(nd,0);
// 		fPrescribedField(i,1) = scale*0.002*initial_coordinates(nd,0);
		double x = initial_coordinates(nd,0);
		double y = initial_coordinates(nd,1);

 		// Linear displacement field		
//  		double A1 =  0.010;
//  		double B1 =  0.002;
//  		double C1 = -0.000;
//  		double A2 =  0.015;
//  		double B2 = -0.008;
//  		double C2 =  0.000;
//  		fPrescribedField(i,0) = scale * (A1*x + B1*y + C1*x*y);
//  		fPrescribedField(i,1) = scale * (A2*x + B2*y + C2*x*y);

		// Manufactured solution 1
		double A = 1.0;
		double B = 1.0;
		double C = 1.0;
		double D = 1.0;
		double Pi = acos(-1.0);
		double U = sin(A*Pi*x)*cos(B*Pi*y);
		double V = sin(C*Pi*x)*cos(D*Pi*y);
 		fPrescribedField(i,0) = scale*U;
 		fPrescribedField(i,1) = scale*V;
		
	}
	
	/* compute forces based on error */
	int numBCs = fConstrainedDOFs.Length();
 	int valueIndex = 0;
 	fConstraintValues = 0.0;
 	for (int i = 0; i < numBCs; i++) {	
		int dof = fConstrainedDOFs[i];
		int nnd = fNodeSets.MinorDim(i);
		for (int j = 0; j < nnd; j++)
			fConstraintValues[valueIndex++] = -kforce*(us(j, dof) - fPrescribedField(j, dof));
	}
}

#endif /* CONTINUUM_ELEMENT */
