/* $Id: NodeManagerT.cpp,v 1.62.2.5 2005-06-10 23:02:01 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "NodeManagerT.h"
#include "ElementsConfig.h"

#include <iostream.h>
#include <iomanip.h>
#include <limits.h>
#include <ctype.h>

#include "ifstreamT.h"
#include "FEManagerT.h"
#include "IOManager.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "LocalArrayT.h"
#include "nIntegratorT.h"
#include "eIntegratorT.h"
#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "PartitionT.h"
#include "ReLabellerT.h"
#include "OutputSetT.h"
#include "ParameterUtils.h"
#include "FieldT.h"

/* force BC controllers */
#include "AugLagWallT.h"
#include "PenaltyWallT.h"
#include "PenaltySphereT.h"
#include "AugLagSphereT.h"
#include "MFPenaltySphereT.h"
#include "PenaltyCylinderT.h"
#include "AugLagCylinderT.h"
#include "MFAugLagMultT.h"
#include "FieldMFAugLagMultT.h"

/* kinematic BC controllers */
#include "K_FieldT.h"
#include "BimaterialK_FieldT.h"
#include "MappedPeriodicT.h"
#include "TiedNodesT.h"
#include "PeriodicNodesT.h"
#include "ScaledVelocityNodesT.h"
#include "SetOfNodesKBCT.h"
#include "TorsionKBCT.h"
#include "ConveyorT.h"

using namespace Tahoe;

/* constructor */
NodeManagerT::NodeManagerT(FEManagerT& fe_manager, CommManagerT& comm_manager):
	ParameterInterfaceT("nodes"),
	fFEManager(fe_manager),
	fCommManager(comm_manager),
	fCoordUpdateIndex(-1),
	fInitCoords(NULL),
	fCoordUpdate(NULL),
	fCurrentCoords(NULL),
	fNeedCurrentCoords(false)
{
	/* set console */
	iSetName("nodes");
	
	/* init support */
	fFieldSupport.SetFEManager(&fe_manager);
	fFieldSupport.SetNodeManager(this);	
}

/* destructor */
NodeManagerT::~NodeManagerT(void) 
{ 
	/* free fields */
	for (int i = 0; i < fFields.Length(); i++)
		delete fFields[i];

	/* free current coordinate array */
	if (fCurrentCoords != fInitCoords) delete fCurrentCoords;
}

/* basic MP support */
int NodeManagerT::Rank(void) const { return fFEManager.Rank(); }
int NodeManagerT::Size(void) const { return fFEManager.Size(); }

int NodeManagerT::NumEquations(int group) const 
{
	/* all fields store the same number of equations */
	int neq = 0;
	for (int i = 0; neq == 0 && i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			neq += fFields[i]->NumEquations();
			
	/* just XDOF equations? */
	if (neq == 0) neq += XDOF_ManagerT::NumEquations(group);
				
	return neq; 
}

int NodeManagerT::NumFields(int group) const 
{
	int num_fields = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			num_fields++;
	return num_fields; 
}

int NodeManagerT::NumDOF(int group) const
{
	/* sum over fields */
	int ndof = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			ndof += fFields[i]->NumDOF();
			
	return ndof; 
}

/* return a pointer to the specified load time function */
const ScheduleT* NodeManagerT::Schedule(int num) const
{
	return fFEManager.Schedule(num);
}

/* register data for output */
void NodeManagerT::RegisterOutput(void)
{
	/* register output from fields */
	for (int j = 0; j < fFields.Length(); j++)
		fFields[j]->RegisterOutput();

	/* configure output for history node sets */
	int num_sets = fHistoryNodeSetIDs.Length();
	if (num_sets > 0)
	{
		/* history nodes written on per field for every set */
		int num_ID = num_sets*fFields.Length();
	    fHistoryOutputID.Dimension(num_ID,3);

		/* configure output of history nodes */
		ModelManagerT& model = *(fFEManager.ModelManager());
		int dex = 0;
		for (int j = 0; j < fFields.Length(); j++) /* loop over fields */
		{
			FieldT& field = *(fFields[j]);
		
			/* field labels */
			const ArrayT<StringT>& labels = field.Labels();
			int ndof = labels.Length();
			
			/* output labels */
			ArrayT<StringT> n_labels(((field.Order() + 1) + 1)*ndof); /* all derivatives + force */

			/* loop over time derivatives */
			int label_dex = 0;
			StringT suffix = "_";
			for (int k = 0; k <= field.Order(); k++)
			{
				for (int l = 0; l < ndof; l++)
				{
					n_labels[label_dex] = labels[l];
					if (k > 0) n_labels[label_dex].Append(suffix);
					label_dex++;
				}
				suffix.Append("t");
			}

			for (int i = 0; i < ndof; i++) /* force */
				n_labels[label_dex++].Append("F_", n_labels[i]);

			for (int i = 0; i < num_sets; i++) /* loop over id sets */
			{
				/* set identifier */
				const StringT& ID = fHistoryNodeSetIDs[i];
		
				/* specify output - "free set" */
				OutputSetT output_set(model.ElementGroupGeometry(ID), model.ElementGroup(ID), n_labels);
				 
				/* register output */
				fHistoryOutputID(dex,0) = fFEManager.RegisterOutput(output_set);
				fHistoryOutputID(dex,1) = j;
				fHistoryOutputID(dex,2) = i;
			
				/* register the node set as a "connectivity" */
				if (j == 0) /* once for each set */
				{
					/* reorder "connectivities" as nodes used - do this because "nodal"
				 	 * output for the set must be written in the order of nodes used and
				 	 * original "connectivities" are not guaranteed to be in this order */
					const iArrayT& nodes_used = output_set.NodesUsed();
					iArray2DT new_conn(nodes_used.Length(), 1, nodes_used.Pointer());
					model.UpdateElementGroup(ID, new_conn, false);
				}
				
				/* next output set */
				dex++;
			}
		}
	}
}

/* collect fields with the given group ID */
void NodeManagerT::CollectFields(int group, ArrayT<FieldT*>& fields) const
{
	/* count */
	int count = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			count++;

	fields.Dimension(count);
	if (count == 0) 
		return;
	else
	{
		/* collect */
		count = 0;
		for (int i = 0; i < fFields.Length(); i++)
			if (fFields[i]->Group() == group)
				fields[count++] = fFields[i];	
	}
}

void NodeManagerT::Equations(int group, AutoArrayT<const iArray2DT*>& eq_1, 
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* from fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->EquationSets(eq_1, eq_2);
}

void NodeManagerT::ConnectsU(int group, 
	AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
	/* from fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->Connectivities(connects_1, connects_2, equivalent_nodes);
}

/* return the implicit-explicit flag for the given group */
IntegratorT::ImpExpFlagT NodeManagerT::ImplicitExplicit(int group) const
{
	IntegratorT::ImpExpFlagT flag = IntegratorT::kExplicit;

	/* loop over fields in the group */
	for (int i = 0; flag == IntegratorT::kExplicit && i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			flag = fFields[i]->nIntegrator().ImplicitExplicit();

	return flag;
}

/* return a pointer to the field with the specified name */
const FieldT* NodeManagerT::Field(const char* name) const
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i] && fFields[i]->FieldName() == name)
			return fFields[i];

	/* not found */
	return NULL;
}

/* register the local coordinate array with its source */
void NodeManagerT::RegisterCoordinates(LocalArrayT& array) const
{
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
			NodeManagerT* non_const_this = (NodeManagerT*) this;
			non_const_this->fNeedCurrentCoords = true;
			break;					
		}
		default:
			ExceptionT::GeneralFail("NodeManagerT::RegisterCoordinates", 
				"not a coordinate type: %d", array.Type());
	}
}

CommManagerT& NodeManagerT::CommManager(void) const { return fCommManager; }

/* read/write access to the coordinate update field */
dArray2DT* NodeManagerT::CoordinateUpdate(void)
{
	if (!fCoordUpdate)
		return NULL;
	else
		return &((*fCoordUpdate)[0]); /* zeroth order component */
}

GlobalT::SystemTypeT NodeManagerT::TangentType(int group) const
{
	/* initialize */
	GlobalT::SystemTypeT type = GlobalT::kUndefined;

	/* check field in this group */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			type = GlobalT::MaxPrecedence(type, fFields[i]->SystemType());
	return type;
}

/* (re-)set the system configuration */
GlobalT::InitStatusT NodeManagerT::UpdateConfiguration(int group)
{
	GlobalT::InitStatusT status = GlobalT::kContinue;
	
	/* equations have never been assigned */
	int neq = NumEquations(group);
	if (neq < 0) status = GlobalT::MaxPrecedence(status, GlobalT::kNewEquations);
//NOTE: should this be kNewEquations or kNewInteractions? kNewEquations will trigger
//      a repeat of the loop. Here we are saying equations need to be re-assigned, but
//      not because we have changed the number of unknowns.

#pragma message("UpdateConfiguration for FBC and KBC controllers?")
	
	return status;
}

/* apply kinematic boundary conditions */
void NodeManagerT::InitStep(int group)
{
	const char caller[] = "NodeManagerT::InitStep";

	/* decomposition information */
	PartitionT::DecompTypeT decomp = fCommManager.DecompType();
	const ArrayT<int>* partition_nodes = fCommManager.PartitionNodes();
	if (decomp == PartitionT::kIndex && !partition_nodes)
		ExceptionT::GeneralFail(caller, "expecting non-NULL partition nodes");

	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
	{
		/* initialize for current group */
		if (fFields[i]->Group() == group)
		{
			/* update bounds */
			int& field_start = fFieldStart[i];
			int& field_end = fFieldEnd[i];
			if (false && decomp == PartitionT::kIndex) {

				/* range of (real) nodes updated by this processor */
				int beg = partition_nodes->First();
				int end = partition_nodes->Last();
	
				/* set limits */
				int ndof = fFields[i]->NumDOF();
				field_start = ndof*beg;
				field_end = (ndof*end) + (ndof - 1);	
			}

			/* initialize step */
			fFields[i]->InitStep(field_start, field_end);
		}
	}

	/* update current configurations */
	if (fCoordUpdate && fCoordUpdate->Group() == group) UpdateCurrentCoordinates();

	/* clear history of relaxation over tbe last step */
	fXDOFRelaxCodes[group] = GlobalT::kNoRelax;
}

/* compute the nodal contribution to the tangent */
void NodeManagerT::FormLHS(int group, GlobalT::SystemTypeT sys_type)
{
	/* skip for explicit dynamics */
	for (int i = 0; i < fFields.Length(); i++) {
		FieldT* field = fFields[i];
		if (field->Group() == group && field->Integrator().ImplicitExplicit() == IntegratorT::kImplicit)
			field->FormLHS(sys_type);
	}
}
	
/* compute the nodal contribution to the residual force vector */
void NodeManagerT::FormRHS(int group)
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->FormRHS();
}

/* call to signal end of RHS calculation to allow NodeManagerT to post-process
 * the total system force */
void NodeManagerT::EndRHS(int group)
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->EndRHS();
}

/* call to signal end of LHS calculation to allow NodeManagerT to post-process
 * the total system tangent matrix */
void NodeManagerT::EndLHS(int group)
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->EndLHS();
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT NodeManagerT::RelaxSystem(int group)
{
	/* initialize */
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;

	/* check fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			relax = GlobalT::MaxPrecedence(relax, fFields[i]->RelaxSystem());

	/* check external DOF groups */
	bool reset_XDOF = XDOF_ManagerT::ResetTags(group);
	GlobalT::RelaxCodeT code_XDOF = (reset_XDOF) ? GlobalT::kReEQ : GlobalT::kNoRelax;

	return GlobalT::MaxPrecedence(relax, code_XDOF);
}

/* set the time step */
void NodeManagerT::SetTimeStep(double dt) {
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->SetTimeStep(dt);
}

/* update the active degrees of freedom */
void NodeManagerT::Update(int group, const dArrayT& update)
{
	/* update fields */
	for (int i = 0; i < fFields.Length(); i++)
	{
		if (fFields[i]->Group() == group)
		{
			/* assemble contribution from local solver */
			fFields[i]->AssembleUpdate(update);

			/* gather/distribute external contribution */
			fCommManager.AllGather(fMessageID[i], fFields[i]->Update());
			
			/* apply the update */
			fFields[i]->ApplyUpdate(fFieldStart[i], fFieldEnd[i]);
		}
	}

	/* update current configurations */
	if (fCoordUpdate && fCoordUpdate->Group() == group)
		UpdateCurrentCoordinates();
	
	/* inherited - update external DOF */
	XDOF_ManagerT::Update(group, update);
}

/* update the current configuration. This is called by NodeManagerT::Update
	 * and does not usually need to be called explicitly. */
void NodeManagerT::UpdateCurrentCoordinates(void)
{
	const char caller[] = "NodeManagerT::UpdateCurrentCoordinates";
	if (fCoordUpdate)
	{
		/* should be allocated */
		if (!fCurrentCoords) ExceptionT::GeneralFail(caller, "current coords not initialized");

		/* bounds */
		int field_start = fFieldStart[fCoordUpdateIndex];
		int field_end = fFieldEnd[fCoordUpdateIndex];
	
		/* simple update assuming displacement degrees of freedom are the
		 * nodal values */
		if (field_end >= field_start)
		{
			/* update coordinates of owned nodes */
			fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0], field_start, field_end);
			
			/* update coordinates of image nodes */
			int first_ghost = fCommManager.NumRealNodes();
			int offset = first_ghost*NumSD();
			double* px = fCurrentCoords->Pointer(offset);
			const double* pX = InitialCoordinates().Pointer(offset);
			const double* pu = (*fCoordUpdate)[0].Pointer(offset);
			int len = fCurrentCoords->Length();
			for (int i = offset; i < len; i++)
				*px++ = (*pX++) + (*pu++);
			/* NOTE: could use CommManagerT::GhostNodes and CommManagerT::NodesWithGhosts
			 *       to compute current coordinates only for images of nodes owned by this
			 *       processor and then exchange to update all other processors. However,
			 *       CommManagerT::AllGather assumes an exchange of the full field (all nodes)
			 *       and then copies values to image nodes, so another strict (no automatic
			 *       handling of ghost node values), decomposition-independent exchange
			 *       method would neeed to be added to CommManagerT.
			 */
		}
		else if (field_end == -1)
			fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
		else
			ExceptionT::GeneralFail(caller, "field_end has unexpected value %d", field_end);
	}	
}

/* update history */
void NodeManagerT::CloseStep(int group)
{
	/* loop over fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->CloseStep();
}

/* initial condition/restart functions */
void NodeManagerT::InitialCondition(void)
{
	/* relaxation flags */
	fXDOFRelaxCodes.Dimension(fFEManager.NumGroups());
	fXDOFRelaxCodes = GlobalT::kNoRelax;

	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
	{
		FieldT& field = *(fFields[i]);

		/* apply initial conditions */
		field.InitialCondition();

		/* gather/distribute external contribution */
		for (int j = 0; j <= field.Order(); j++)
			fCommManager.AllGather(fMessageID[i], field[j]);
	}

	/* update current configurations */
	UpdateCurrentCoordinates();
}

void NodeManagerT::ReadRestart(ifstreamT& in)
{
	/* nodes owned by this partition */
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();

	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
	{
		FieldT& field = *(fFields[i]);
		field.ReadRestart(in, part_nodes);
		
		/* gather/distribute external contribution */
		for (int j = 0; j <= field.Order(); j++)
			fCommManager.AllGather(fMessageID[i], field[j]);
		
		/* reset history */
		field.CloseStep();
	}
	
	/* update current configurations */
	UpdateCurrentCoordinates();
}

void NodeManagerT::WriteRestart(ofstreamT& out) const
{
	/* nodes owned by this partition */
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();

	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->WriteRestart(out, part_nodes);
}

/* reset displacements (and configuration to the last known solution) */
GlobalT::RelaxCodeT NodeManagerT::ResetStep(int group)
{
	/* initialize return value */
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;

	/* reset fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			relax = GlobalT::MaxPrecedence(relax, fFields[i]->ResetStep());

	/* reset the XDOF elements */
	XDOF_ManagerT::ResetState(group);

	/* inherited - reset external DOF */
	return GlobalT::MaxPrecedence(relax, XDOF_ManagerT::ResetTags(group));
}

void NodeManagerT::WriteOutput(void)
{
	/* stream */
	ostream& out = fFEManager.Output();

	/* loop over fields */
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->WriteOutput(out);

	/* external DOF */
	if (fDOFElements.Length() > 0)
	{
		ostream& out = fFEManager.Output();
	
		out << "\n E l e m e n t   d e g r e e s   o f   f r e e d o m :\n\n";
		out << " Number of element equation groups . . . . . . . = ";
		out << fDOFElements.Length() << "\n\n";
	
		int set_index = -1;
		for (int i = 0 ; i < fDOFElements.Length(); i++)
		{
			out << " Group " << i+1 << ":\n";
			for (int j = 0; j < fNumTagSets[i]; j++)
			{
				out << " Set " << j+1 << ":\n";
				set_index++;
				WriteData(out, "Element degrees of freedom", "dof",
					*(fXDOFs[set_index]), &(fDOFElements[i]->DOFTags(j)));
			}
		}
	}

#pragma message("NodeManagerT -- Necessary kludge here")
	UpdateCurrentCoordinates();

	/* nodal histories */
//NOTE - moved here from CloseStep	
	WriteNodalHistory();
}

void NodeManagerT::SetEquationNumbers(int group)
{
	const char caller[] = "NodeManagerT::SetEquationNumbers";
	int token = 1;
	try {
	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0)
		ExceptionT::GeneralFail(caller, "group has no fields: %d", group);
	
	/* initialize equations numbers arrays */
	for (int i = 0; i < fields.Length(); i++)
		fields[i]->InitEquations();

	/* mark external nodes */
	const ArrayT<int>* ex_nodes = fCommManager.ExternalNodes();
	if (ex_nodes)
		for (int i = 0; i < fields.Length(); i++)
		{
			/* field equations array */
			iArray2DT&  eqnos = fields[i]->Equations();

			/* mark all external as inactive for setting local
			 * equation numbers */
			for (int j = 0; j < ex_nodes->Length(); j++)
				eqnos.SetRow((*ex_nodes)[j], FieldT::kExternal);	
		}

	/* assign active equation numbers node-by-node across fields
	 * in the group */
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
	int num_eq = 0;
	int nnd = (part_nodes) ? part_nodes->Length() : NumNodes();
	for (int i = 0; i < nnd; i++)
		for (int j = 0; j < fields.Length(); j++)
		{
			int ndof = fields[j]->NumDOF();
			int nd = (part_nodes) ? (*part_nodes)[i] : i;
			int* peq = (fields[j]->Equations())(nd);
			for (int k = 0; k < ndof; k++)
			{
				/* active equation */
				if (*peq >= FieldT::kInit) 
					*peq = ++num_eq;
	
				peq++;
			}
		}

	/* assign equation numbers to XDOF's */
	XDOF_ManagerT::SetEquations(group, num_eq);

	/* set equation arrays - assume all start at 1 */
	int start_eq = 1;
	for (int i = 0; i < fields.Length(); i++)
		fields[i]->FinalizeEquations(start_eq, num_eq);
	} /* end try */
	
	catch (ExceptionT::CodeT error) { token = 0; }

	/* synch and check */
	if (fCommManager.Communicator().Sum(token) != Size()) ExceptionT::BadHeartBeat(caller);
}

void NodeManagerT::RenumberEquations(int group, 
	const ArrayT<const iArray2DT*>& connects_1,
	const ArrayT<const RaggedArray2DT<int>*>& connects_2)
{
	cout << "\n NodeManagerT::RenumberEquations: start" << endl;

	/* bandwidth reducer */
	ReLabellerT relabel;

	/* send to relabeller */
	for (int j = 0; j < connects_1.Length(); j++)
		relabel.AddGroup(*(connects_1[j]));
	for (int k = 0; k < connects_2.Length(); k++)
		relabel.AddGroup(*(connects_2[k]));
	
	/* collect sets of equation numbers */
	AutoArrayT<iArray2DT*> eqnos;
	EquationNumbers(group, eqnos);

	int numtest = relabel.Renumber(eqnos);
	if (numtest != NumEquations(group))
		ExceptionT::GeneralFail("NodeManagerT::RenumberEquations", 
			"expecting to renumber %d eqns, but hit %d", NumEquations(group), numtest);
	
	/* rearrange equations if needed */
	CheckEquationNumbers(group);

	/* reset fields */
	for (int j = 0; j < fFields.Length(); j++)
		if (fFields[j]->Group() == group)
		{
			int start = fFields[j]->EquationStart();
			int num_eq = NumEquations(group);
			fFields[j]->FinalizeEquations(start, num_eq);
		}

	cout << "\n NodeManagerT::RenumberEquations: done" << endl;
}

void NodeManagerT::SetEquationNumberScope(int group, GlobalT::EquationNumberScopeT scope)
{
	/* id's of external nodes */
	const ArrayT<int>* ex_nodes = fCommManager.ExternalNodes();

	//TEMP - external DOF's no tested with other scopes
	if (scope != GlobalT::kLocal && ex_nodes && NumTagSets() > 0)
		ExceptionT::GeneralFail("NodeManagerT::SetEquationNumberScope", 
			"external DOF only verified with local numbering");

	/* switch numbering scope - with external nodes */
	if (scope == GlobalT::kGlobal && ex_nodes)
	{
		/* shift local equation numbers */
		int start = fFEManager.GlobalEquationStart(group);
		int shift = start - 1;
		
		/* change numbering scope */
		for (int j = 0; j < fFields.Length(); j++)
			if (fFields[j]->Group() == group)
			{
				iArray2DT& eqnos = fFields[j]->Equations();
				int* peq = eqnos.Pointer();
				int  len = eqnos.Length();
				for (int i = 0; i < len; i++)
				{
					if (*peq > FieldT::kInit) *peq += shift;
					peq++;
				}
				
				/* set up exchange */
				int id = fCommManager.Init_AllGather(eqnos);

				/* gather external contribution */
				fCommManager.AllGather(id, eqnos);

				/* clear exchange */
				fCommManager.Clear_AllGather(id);
			}
			
		/* reset fields */
		for (int j = 0; j < fFields.Length(); j++)
			if (fFields[j]->Group() == group)
			{
				int num_eq = fFields[j]->NumEquations();
				fFields[j]->FinalizeEquations(start, num_eq);
			}
	}
}

/* equation for the given degree of freedom. Equations > 0 are "active" */
int NodeManagerT::EquationNumber(int field, int node, int dof) const
{
	const iArray2DT& eqnos = fFields[field]->Equations();
	return eqnos(node, dof);
}

void NodeManagerT::WriteEquationNumbers(int group, ostream& out) const
{
	/* collect fields in this group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);

	/* print header */
	out << "\n N o d a l   E q u a t i o n   N u m b e r s :\n\n";
	out << " Number of element equation groups . . . . . . . = " << fFEManager.NumGroups() << '\n';
	out << " Group number. . . . . . . . . . . . . . . . . . = " << group+1 << '\n';	
	out << " Number of fields. . . . . . . . . . . . . . . . = " << fields.Length() << '\n';
	
	/* equations per field */
	for (int i = 0; i < fields.Length(); i++)
		fields[i]->WriteEquationNumbers(out, fCommManager.NodeMap());
	
	/* external equation groups */
	for (int i = 0; i < fXDOF_Eqnos.Length(); i++)
		if (fDOFElements[i] -> Group() == group)
		{
			out << "\n XDOF equation set: " << i+1 << '\n';
			fXDOF_Eqnos[i]->WriteNumbered(out);
		}	
}

/* return the current values of the unknowns */
void NodeManagerT::GetUnknowns(int group, int order, dArrayT& unknowns) const
{
	const char caller[] = "NodeManagerT::GetUnknowns";

//TEMP - not fully implemented
if (NumTagSets() > 0) ExceptionT::GeneralFail(caller, "not implemented for XDOF tags");

	/* check */
	int num_eq = NumEquations(group);
	if (unknowns.Length() != num_eq) ExceptionT::GeneralFail(caller);

	/* loop over groups */
	int checksum = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
		{
			/* field data */
			FieldT& field = *(fFields[i]);
			
			/* check order of field */
			if (field.Order() < order)
				ExceptionT::OutOfRange(caller, "order %d is out of range {0,%d}", order, field.Order());

			const dArray2DT& u = field[order];
			const iArray2DT& eqnos = field.Equations();
		
			/* fill values from field */
			const int*   peq = eqnos.Pointer();
			const double* pu = u.Pointer();
			int len = eqnos.Length();
			for (int j = 0; j < len; j++)
			{
				int eq = *peq++;
				if (eq > FieldT::kInit)
				{
					unknowns[eq - 1] = *pu;
					checksum++;
				}
				pu++;
			}
		}

	/* check */
	if (checksum != num_eq) ExceptionT::GeneralFail(caller, "checksum error");
}

/* weight the computational effort of every node */
void NodeManagerT::WeightNodalCost(iArrayT& weight) const
{
	int* p = weight.Pointer();
	int  n = weight.Length();
	for (int i = 0; i < n; i++)
	{
		if (*p < 1) *p = 1;
		p++;
	}	
}

/* reset the number of nodes */
void NodeManagerT::ResizeNodes(int num_nodes)
{
	/* reference coordinates */
	fFEManager.ModelManager()->ResizeNodes(num_nodes);
#pragma message("resize reference coords here or require separate call?")

	/* current coordinates */
	if (fCurrentCoords) fCurrentCoords_man.SetMajorDimension(num_nodes, true);

	/* resize fields */
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->Dimension(num_nodes, true);

	/* averaging work space */
	SetNumAverageRows(NumNodes());
}

/* copy nodal information */
void NodeManagerT::CopyNodeToNode(const ArrayT<int>& source, 
	const ArrayT<int>& target)
{
	/* check */
	if (source.Length() != target.Length()) 
		ExceptionT::SizeMismatch("NodeManagerT::CopyNodeToNode");

	/* copy fields */
	for (int i = 0; i < fFields.Length(); i++)
	{
		/* the field */
		FieldT& field = *(fFields[i]);

		/* copy field data */
		field.CopyNodeToNode(source, target);

		/* gather/distribute external contribution */
		for (int j = 0; j <= field.Order(); j++)
			fCommManager.AllGather(fMessageID[i], field[j]);
		
		/* reset history */
		field.CloseStep();
	}

	/* reset current coordinates */
	UpdateCurrentCoordinates();
}

/* size of the nodal package */
int NodeManagerT::PackSize(void) const
{
	int size = 0;
	for (int i = 0; i < fFields.Length(); i++) {
		FieldT& field = *(fFields[i]);
		size += 2*field.NumDOF()*(field.Order() + 1);
	}
	return size;
}

/* copy field information into the array */
void NodeManagerT::Pack(int node, dArrayT& values) const
{
	int index = 0;
	for (int i = 0; i < fFields.Length(); i++) /* loop over fields */
	{
		FieldT& field = *(fFields[i]);
		int order = field.Order();
		int ndof  = field.NumDOF();
		if (values.Length() < index + 2*ndof*(order + 1))
			ExceptionT::SizeMismatch("NodeManagerT::Pack");
	
			/* loop over time derivatives */
			for (int i = 0; i < order+1; i++) 
			{
				/* values from last step */
				dArray2DT& f = field(0,i);
				dArray2DT& f_last = field(-1,i); 

				/* values at current time */
				f.RowCopy(node, values.Pointer(index));
				index += ndof;

				/* values from the last time step */
				f_last.RowCopy(node, values.Pointer(index));
				index += ndof;
			}	
	}
}

/* write information from the array into the fields */
void NodeManagerT::Unpack(int node, dArrayT& values)
{
	int index = 0;
	for (int i = 0; i < fFields.Length(); i++) /* loop over fields */
	{
		FieldT& field = *(fFields[i]);
		int order = field.Order();
		int ndof  = field.NumDOF();
		if (values.Length() < index + 2*ndof*(order + 1))
			ExceptionT::SizeMismatch("NodeManagerT::Unpack");
	
			/* loop over time derivatives */
			for (int i = 0; i < order+1; i++) 
			{
				/* values from last step */
				dArray2DT& f = field(0,i);
				dArray2DT& f_last = field(-1,i);

				/* values at current time */
				f.SetRow(node, values.Pointer(index));
				index += ndof;

				/* values from the last time step */
				f_last.SetRow(node, values.Pointer(index));
				index += ndof;
			}	
	}
}

//TEMP - trap parallel execution with XDOF
void NodeManagerT::XDOF_Register(DOFElementT* group, const iArrayT& numDOF)
{
	//TEMP - parallel execution not yet supported
	if (fFEManager.Size() > 1)
		ExceptionT::GeneralFail("NodeManagerT::XDOF_Register", "not for parallel execution");
//NOTE: to parallelize XDOF:
// (1) analyze external nodes to see if they interact with any element-generated DOF's
// (2) collect and send these tags/equation numbers separate from primary variables	

	/* inherited */
	XDOF_ManagerT::XDOF_Register(group, numDOF);
}

/* collection equation numbers for mixed connectivities. */
void NodeManagerT::XDOF_SetLocalEqnos(int group, const iArrayT& nodes, 
	iArray2DT& eqnos)
{
	const char caller[] = "NodeManagerT::XDOF_SetLocalEqnos";

	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0)
		ExceptionT::GeneralFail(caller, "group %d has no fields", group);

	/* dimensions */
	int nnd = NumNodes();
	int nen = nodes.Length();
	int neq = eqnos.Length();

	const int* ien = nodes.Pointer();
	int* peq = eqnos.Pointer();
		
	/* count assigned equation numbers */
	int eq_count = 0;
	
	/* loop over element tags */
	for (int j = 0; j < nen; j++)
	{
		int tag = *ien++;
		int tag_offset = 0;
		
		/* loop over fields if needed */
		int dex = 0;
		bool done = false;
		while (!done)
		{
			/* source for equations */
			const iArray2DT* eqnos_source;
		
			/* node tag */
			if (tag < nnd) 
			{
				eqnos_source = &(fields[dex++]->Equations());
				done = (dex == fields.Length()); /* loop over fields */
			}
			else /* XDOF tag */
			{
				/* no loop over fields */
				done = true;
				
				/* resolve tag into its set */
				int tag_set;
				if (!ResolveTagSet(tag, tag_set, tag_offset))
					ExceptionT::GeneralFail(caller, "could not resolve tag into set %d", tag);

				/* equations from tag set */
				eqnos_source = fXDOF_Eqnos[tag_set];
			}		
			
			/* dimension */
			int ndof = eqnos_source->MinorDim();
			eq_count += ndof;

			/* check number of assigned equations */
			if (eq_count > neq)
				ExceptionT::SizeMismatch(caller, "error assigning equations");
		
			/* copy equations */
			eqnos_source->RowCopy(tag - tag_offset, peq);
			
			/* next */
			peq += ndof;			
		}
	}
}

/* collection equation numbers for mixed connectivities. */
void NodeManagerT::XDOF_SetLocalEqnos(int group, const iArray2DT& nodes, 
	iArray2DT& eqnos) const
{
	const char caller[] = "NodeManagerT::XDOF_SetLocalEqnos";

	/* check */
	if (nodes.MajorDim() != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);

	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0)
		ExceptionT::GeneralFail(caller, "group %d has no fields", group+1);

	/* dimensions */
	int nnd = NumNodes();
	int nel = nodes.MajorDim();
	int nen = nodes.MinorDim();
	int neq = eqnos.MinorDim();

	const int* ien = nodes.Pointer();
	int* peq = eqnos.Pointer();
	for (int i = 0; i < nel; i++)
	{
		/* count assigned equation numbers */
		int eq_count = 0;
	
		/* loop over element tags */
		for (int j = 0; j < nen; j++)
		{
			int tag = *ien++;
			int tag_offset = 0;
	
			/* loop over fields if needed */
			int dex = 0;
			bool done = false;
			while (!done)
			{
				/* source for equations */
				const iArray2DT* eqnos_source;

				/* node tag */
				if (tag < nnd)
				{ 
					eqnos_source = &(fields[dex++]->Equations());
					done = (dex == fields.Length()); /* loop over fields */
				}
				else /* XDOF tag */
				{
					/* no loop over fields */
					done = true;
				
					/* resolve tag into its set */
					int tag_set;
					if (!ResolveTagSet(tag, tag_set, tag_offset))
						ExceptionT::GeneralFail(caller, "could not resolve tag %d into set", tag);

					/* equations from tag set */
					eqnos_source = fXDOF_Eqnos[tag_set];
				}

				/* dimension */
				int ndof = eqnos_source->MinorDim();
				eq_count += ndof;

				/* check number of assigned equations */
				if (eq_count > neq)
					ExceptionT::SizeMismatch(caller, "error assigning equations");

				/* copy equations */
				eqnos_source->RowCopy(tag - tag_offset, peq);
			
				/* next */
				peq += ndof;
			}
		}
		
		/* check */
		if (eq_count != neq)
			ExceptionT::GeneralFail(caller, "expecting %d equations not %a", neq, eq_count);
	}
}

/* collection equation numbers for mixed connectivities. */
void NodeManagerT::XDOF_SetLocalEqnos(int group, const RaggedArray2DT<int>& nodes, 
	RaggedArray2DT<int>& eqnos) const
{
	const char caller[] = "NodeManagerT::XDOF_SetLocalEqnos";

	/* check */
	if (nodes.MajorDim() != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);

	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0)
		ExceptionT::GeneralFail(caller, "group %d has no fields", group+1);

	/* dimensions */
	int nnd = NumNodes();
	int nel = nodes.MajorDim();

	const int* ien = nodes.Pointer();
	int* peq = eqnos.Pointer();
	for (int i = 0; i < nel; i++)
	{
		/* dimensions */
		int nen = nodes.MinorDim(i);
		int neq = eqnos.MinorDim(i);

		/* count assigned equation numbers */
		int eq_count = 0;
	
		/* loop over element tags */
		for (int j = 0; j < nen; j++)
		{
			int tag = *ien++;
			int tag_offset = 0;
	
			/* loop over fields if needed */
			int dex = 0;
			bool done = false;
			while (!done)
			{
				/* source for equations */
				const iArray2DT* eqnos_source;

				/* node tag */
				if (tag < nnd)
				{ 
					eqnos_source = &(fields[dex++]->Equations());
					done = (dex == fields.Length()); /* loop over fields */
				}
				else /* XDOF tag */
				{
					/* no loop over fields */
					done = true;
				
					/* resolve tag into its set */
					int tag_set;
					if (!ResolveTagSet(tag, tag_set, tag_offset))
						ExceptionT::GeneralFail(caller, "could not resolve tag %d into set", tag);

					/* equations from tag set */
					eqnos_source = fXDOF_Eqnos[tag_set];
				}

				/* dimension */
				int ndof = eqnos_source->MinorDim();
				eq_count += ndof;

				/* check number of assigned equations */
				if (eq_count > neq)
					ExceptionT::SizeMismatch(caller, "error assigning equations");

				/* copy equations */
				eqnos_source->RowCopy(tag - tag_offset, peq);
				
				/* next */
				peq += ndof;
			}
		}

		/* check */
		if (eq_count != neq)
			ExceptionT::GeneralFail(caller, "expecting %d equations not %d", neq, eq_count);
	}
}

/* describe the parameters needed by the interface*/
void NodeManagerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* name of the field which updates the coordinates */
	ParameterT coord_update(ParameterT::Word, "coordinate_update_field");
	coord_update.SetDefault("displacement");
	list.AddParameter(coord_update, ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameter lists */
void NodeManagerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the fields */
	sub_list.AddSub("field", ParameterListT::OnePlus);
	
	/* list of history node ID's */
	sub_list.AddSub("history_node_ID_list", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NodeManagerT::NewSub(const StringT& name) const
{
	if (name == "field")
		return new FieldT(fFieldSupport);
	else
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void NodeManagerT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "NodeManagerT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* read coordinates information */
	SetCoordinates();

	/* relaxation flags */
	fXDOFRelaxCodes.Dimension(fFEManager.NumGroups());
	fXDOFRelaxCodes = GlobalT::kNoRelax;

	/* collect history nodes */
	const ParameterListT* history_nodes = list.List("history_node_ID_list");
	if (history_nodes) {
	
		/* model database */
		ModelManagerT* model = fFEManager.ModelManager();
	
		/* external nodes - no output written */
		const ArrayT<int>* p_ex_nodes = fCommManager.ExternalNodes();
		iArrayT ex_nodes;
		if (p_ex_nodes) ex_nodes.Alias(*p_ex_nodes);

		int num_ID = history_nodes->NumLists("String");
		fHistoryNodeSetIDs.Dimension(num_ID);
		for (int i = 0; i < num_ID; i++) {
		
			/* read node set */
			const StringT& node_set_ID = history_nodes->GetList("String", i).GetParameter("value");		
			const iArrayT& node_set = model->NodeSet(node_set_ID);
		
			/* slow-but-steady way */
			ArrayT<bool> is_external(node_set.Length());
			is_external = true;
			int count = 0;
			for (int j = 0; j < node_set.Length(); j++)
				if (!ex_nodes.HasValue(node_set[j])) {
					is_external[j] = false;
					count++;
				}
		
			/* collect non-external nodes */
			iArray2DT set(count, 1);
			count = 0;
			for (int j = 0; j < node_set.Length(); j++)
				if (!is_external[j])
					set[count++] = node_set[j];

			/* register the set with the model manager */
			StringT ID = "55";
			ID = model->FreeElementID(ID);
			if (!model->RegisterElementGroup(ID, set, GeometryT::kPoint, true))
				ExceptionT::BadInputValue(caller, "error initializing node set %d as model element ID %d", 
					node_set_ID.Pointer(), ID.Pointer());

			/* store generated ID */
			fHistoryNodeSetIDs[i] = ID;
		}
	}

	/* name of the field which updates the coordinates */
	StringT coords_update_field = "displacement"; /* default update field name - backward compatibility */
	const ParameterT* coord_update = list.Parameter("coordinate_update_field");
	if (coord_update)
		coords_update_field = *coord_update;

	/* construct fields */
	int num_fields = list.NumLists("field");
	fFields.Dimension(num_fields);
	fFields = NULL;
	fMessageID.Dimension(num_fields);
	for (int i = 0; i < fFields.Length(); i++)
	{
		/* parameters */
		const ParameterListT* field_params = list.List("field", i);
	
		/* new field */			
		FieldT* field = new FieldT(fFieldSupport);

		/* store */
		fFields[i] = field;

		/* initialize */
		field->TakeParameterList(*field_params);
		field->Dimension(NumNodes(), false);
		field->Clear();

		/* coordinate update field */
		if (field->FieldName() == coords_update_field) {
			fCoordUpdateIndex = i;
			if (fCoordUpdate) ExceptionT::BadInputValue(caller, "coordinate update field already set");
			fCoordUpdate = field;
			fCurrentCoords = new dArray2DT;
			fCurrentCoords_man.SetWard(0, *fCurrentCoords, NumSD());
			fCurrentCoords_man.SetMajorDimension(NumNodes(), false);
			(*fCurrentCoords) = InitialCoordinates();
		}
			
		/* set up communication of field */
		fMessageID[i] = fCommManager.Init_AllGather(fFields[i]->Update());
	}
	fFieldStart.Dimension(num_fields);
	fFieldEnd.Dimension(num_fields);
	fFieldStart = 0;
	fFieldEnd = -1;
}

/**********************************************************************
 * Protected
 **********************************************************************/

void NodeManagerT::SetCoordinates(void)
{
	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* read coordinates */
	model->ReadCoordinates();
		
	/* set pointer */	
	fInitCoords = &(model->Coordinates());

	/* check element groups to see if node data should be 
	   adjusted to be 2D, some element groups require
	   fNumSD == fNumDOF */
	if (NumSD() == 3 && model->AreElements2D()) {
		cout << "\n NodeManagerT::EchoCoordinates: WARNING: Adjusting nodal data to 2D" << endl;
		model->AdjustCoordinatesto2D();
	}
	  
	/* verbose output */
	if (fFEManager.PrintInput())
	{
		ofstreamT& out = fFEManager.Output();
		int d_width = out.precision() + kDoubleExtra;

		/* print main header */
		out << " Number of nodal points. . . . . . . . . . . . . = " << NumNodes() << '\n';
		out << " Number of spatial dimensions. . . . . . . . . . = " << NumSD() << '\n';
	
		/* write header */
		out << setw(kIntWidth) << "node"
		    << setw(kIntWidth) << "gl.node"
		    << setw(kIntWidth) << "proc";		
		for (int i = 0; i < NumSD(); i++)
			out << setw(d_width - 2) << "x[" << i + 1 << "]";
		out << '\n';

		/* arrays */
		const ArrayT<int>* processor = fCommManager.ProcessorMap();
		const dArray2DT& init_coords = InitialCoordinates();
		const ArrayT<int>* node_map = fCommManager.NodeMap();
		for (int i = 0; i < init_coords.MajorDim(); i++)
		{
			out << setw(kIntWidth) << i+1
			    << setw(kIntWidth) << ((node_map) ? (*node_map)[i]+1 : i+1)
			    << setw(kIntWidth) << ((processor) ? (*processor)[i] : 0);		
			for (int j = 0; j < NumSD(); j++)
				out << setw(d_width) << init_coords(i,j);
			out << '\n';
		}
		out.flush();
	}

	/* set start tag for external DOF */
	XDOF_ManagerT::SetStartTag(NumNodes());
	
	/* averaging work space */
	SetNumAverageRows(NumNodes());
}

/* simple output function */
void NodeManagerT::WriteData(ostream& out, const char* title,
	const char* name, const dArray2DT& data, const iArrayT* rowlabels) const
{
#pragma unused(name)

	int d_width = out.precision() + kDoubleExtra;

	/* data dimension info */
	out << "\n " << title << " :\n\n";	
	out << " Number of nodal points. . . . . . . . . . . . . = " << data.MajorDim() << '\n';
	out << " Number of nodal degrees of freedom. . . . . . . = " << data.MinorDim() << "\n\n";
	
	/* data header */
	out << setw(kIntWidth) << "node";
	for (int i = 1; i <= data.MinorDim(); i++)
		out << setw(d_width - 2) << "d[" << i << "]";
	out << '\n';
	
	/* the data */
	if (rowlabels)
	{
		/* check */
		if (rowlabels->Length() != data.MajorDim()) ExceptionT::SizeMismatch();
	
		for (int i = 0; i < data.MajorDim(); i++)
		{
			out << setw(kIntWidth) << (*rowlabels)[i];
			data.PrintRow(i, out);
		}
	}
	else
		data.WriteNumbered(out);
}

KBC_ControllerT* NodeManagerT::NewKBC_Controller(FieldT& field, int code)
{
	const char caller[] = "NodeManagerT::NewKBC_Controller";
	switch(code)
	{
		case KBC_ControllerT::kPrescribed:
			return new KBC_ControllerT(fFieldSupport);
	
#ifdef CONTINUUM_ELEMENT
		case KBC_ControllerT::kK_Field:
			return new K_FieldT(fFieldSupport);

		case KBC_ControllerT::kBimaterialK_Field:	
			return new BimaterialK_FieldT(fFieldSupport);
#endif

		case KBC_ControllerT::kMappedPeriodic:	
			return new MappedPeriodicT(fFieldSupport, field);

		case KBC_ControllerT::kTiedNodes:
		{
			TiedNodesT* kbc = new TiedNodesT(fFieldSupport, field);
			return kbc;
		}
		case KBC_ControllerT::kPeriodicNodes:
		{
			PeriodicNodesT* kbc = new PeriodicNodesT(fFieldSupport, field);
			return kbc;
		}
		case KBC_ControllerT::kScaledVelocityNodes:
		{
			ScaledVelocityNodesT* kbc = new ScaledVelocityNodesT(fFieldSupport, field);
			return kbc;
		}
		case KBC_ControllerT::kSetOfNodesKBC:
		{
			SetOfNodesKBCT* kbc = new SetOfNodesKBCT(fFieldSupport, field);
			return kbc;
		}
		case KBC_ControllerT::kTorsion:
		{
			TorsionKBCT* kbc = new TorsionKBCT(fFieldSupport);
			return kbc;
		}
#if defined(CONTINUUM_ELEMENT) && defined(COHESIVE_SURFACE_ELEMENT)
		case KBC_ControllerT::kConveyor:
		{
			ConveyorT* kbc = new ConveyorT(fFieldSupport, field);
			return kbc;
		}
#endif

#if 0
		case KBC_ControllerT::kConveyorSym:
		{
			ConveyorSymT* kbc = new ConveyorSymT(fFieldSupport, field);
			return kbc;
		}
#endif

		default:
			ExceptionT::BadInputValue(caller, "KBC controller code %d is not supported", code);
	}
	return NULL;
}

FBC_ControllerT* NodeManagerT::NewFBC_Controller(int code)
{
  const char caller[] = "NodeManagerT::NewFBC_Controller";

	FBC_ControllerT* fbc = NULL;
	switch(code)
	{
		case FBC_ControllerT::kPenaltyWall:
			fbc = new PenaltyWallT;
			break;

		case FBC_ControllerT::kAugLagWall:
			fbc = new AugLagWallT;
			break;

		case FBC_ControllerT::kPenaltySphere:	
			fbc = new PenaltySphereT;
			break;

		case FBC_ControllerT::kPenaltyCylinder:	
			fbc = new PenaltyCylinderT;
			break;

		case FBC_ControllerT::kAugLagSphere:	
			fbc = new AugLagSphereT;
			break;

		case FBC_ControllerT::kMFPenaltySphere:	
			fbc = new MFPenaltySphereT;
			break;

#ifdef CONTINUUM_ELEMENT
	    case FBC_ControllerT::kMFAugLagMult:
	    	fbc = new MFAugLagMultT;
	    	break;

	    case FBC_ControllerT::kFieldMFAugLagMult:
	    	fbc = new FieldMFAugLagMultT;
	    	break;
#endif

	    case FBC_ControllerT::kAugLagCylinder:
	    	fbc = new AugLagCylinderT;
	    	break;

		default:
			ExceptionT::BadInputValue(caller, "FBC controller code %d is not supported", code);
	}	
	return fbc;
}

/* access to global equation numbers */
void NodeManagerT::EquationNumbers(int group, AutoArrayT<iArray2DT*>& equationsets)
{
	/* fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			equationsets.Append(&(fFields[i]->Equations()));

	/* add XDOF equation sets */
	XDOF_ManagerT::EquationNumbers(group, equationsets);
}

void NodeManagerT::CheckEquationNumbers(int group)
{
	/* fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			XDOF_ManagerT::CheckEquationNumbers(fFEManager.Output(), fFields[i]->Equations());
}

/**********************************************************************
 * Private
 **********************************************************************/

void NodeManagerT::WriteNodalHistory(void)
{
	if (fHistoryOutputID.MajorDim() > 0)
	{
		/* model manager */
		ModelManagerT& model = *(fFEManager.ModelManager());

		/* loop over output node sets */
		for (int i = 0; i < fHistoryOutputID.MajorDim(); i++)
		{
			/* set information */
			int ID = fHistoryOutputID(i,0);
			const FieldT& field = *(fFields[fHistoryOutputID(i,1)]);
			int nset = fHistoryOutputID(i,2);
		
			/* node set */
			const iArray2DT& node_set = model.ElementGroup(fHistoryNodeSetIDs[nset]);

			/* conjugate force */
			int ndof = field.NumDOF();
			dArrayT force(ndof);
		
			/* output values */
			dArray2DT n_values(node_set.Length(), ((field.Order() + 1) + 1)*ndof);
			dArray2DT e_values;
		
			/* nodes in set */		
			for (int j = 0; j < node_set.Length(); j++) 
			{
				/* node map resolves processor-local node number to global
				 * node number. No map means serial execution */
				int node = node_set[j];
						
				/* compute reaction force */
				fFEManager.InternalForceOnNode(field, node, force);

				/* loop over time derivatives */
				int dex = 0;
				for (int l = 0; l <= field.Order(); l++)
					for (int k = 0; k < ndof; k++)
						n_values(j, dex++) = field[l](node, k); /* field and derivatives */

				for (int k = 0; k < ndof; k++) /* force */
					n_values(j, dex++) = force[k];
			}
			
			/* send for output */
			fFEManager.WriteOutput(ID, n_values, e_values);
		}
	}
}
