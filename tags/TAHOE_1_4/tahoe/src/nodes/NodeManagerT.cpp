/* $Id: NodeManagerT.cpp,v 1.48 2004-06-17 07:41:49 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "NodeManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <limits.h>
#include <ctype.h>

#include "ifstreamT.h"
#include "FEManagerT.h"
#include "IOManager.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
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
#include "MFAugLagMultT.h"

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
	fFieldSupport(fe_manager, *this),
	fInitCoords(NULL),
	fCoordUpdate(NULL),
	fCurrentCoords(NULL),
	fNeedCurrentCoords(false)
{
	/* set console */
	iSetName("nodes");
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

/* initialize data */
void NodeManagerT::Initialize(void)
{
	/* get streams */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();

	out << "\n N o d a l   D a t a :\n\n";

	/* echo nodal coordinate data */
	EchoCoordinates(in, out);

	/* set fields */
	EchoFields(in, out);

	/* history nodes */
	EchoHistoryNodes(in, out);

	/* relaxation flags */
	fXDOFRelaxCodes.Dimension(fFEManager.NumGroups());
	fXDOFRelaxCodes = GlobalT::kNoRelax;
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
		if (fFields[i] && fFields[i]->Name() == name)
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

/* the local node to home processor map */
const ArrayT<int>* NodeManagerT::ProcessorMap(void) const { return fFEManager.ProcessorMap(); }

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

/* apply kinematic boundary conditions */
void NodeManagerT::InitStep(int group)
{
	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->InitStep();

	/* update current configurations */
	if (fCoordUpdate && fCoordUpdate->Group() == group)
		fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);

	/* clear history of relaxation over tbe last step */
	fXDOFRelaxCodes[group] = GlobalT::kNoRelax;
}

/* compute the nodal contribution to the tangent */
void NodeManagerT::FormLHS(int group, GlobalT::SystemTypeT sys_type)
{
	int analysiscode = fFEManager.Analysis();

//NOTE: this test should really be coming from the controller

	/* skip for explicit dynamics */
	if (analysiscode != GlobalT::kLinExpDynamic   &&
	    analysiscode != GlobalT::kNLExpDynamic    &&
	    analysiscode != GlobalT::kVarNodeNLExpDyn &&
	    analysiscode != GlobalT::kPML)
	{
		for (int i = 0; i < fFields.Length(); i++)
			if (fFields[i]->Group() == group)
				fFields[i]->FormLHS(sys_type);
	}
}
	
/* compute the nodal contribution to the residual force vector */
void NodeManagerT::FormRHS(int group)
{
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->FormRHS();
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

/* update the active degrees of freedom */
void NodeManagerT::Update(int group, const dArrayT& update)
{
	/* update fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
		{
			/* assemble contribution from local solver */
			fFields[i]->AssembleUpdate(update);

			/* gather/distribute external contribution */
			fCommManager.AllGather(fMessageID[i], fFields[i]->Update());
			
			/* apply the update */
			fFields[i]->ApplyUpdate();
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
	if (fCoordUpdate)
	{
		/* should be allocated */
		if (!fCurrentCoords)
			ExceptionT::GeneralFail("NodeManagerT::UpdateCurrentCoordinates", "current coords not initialized");
	
		/* simple update assuming displacement degrees of freedom are the
		 * nodal values */
		fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
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
	if (fCoordUpdate) fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
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
	if (fCoordUpdate) fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
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
	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0)
		ExceptionT::GeneralFail("NodeManagerT::SetEquationNumbers", 
			"group has no fields: %d", group);
	
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
		fields[i]->WriteEquationNumbers(out, fFEManager.NodeMap());
	
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
//TEMP - not fully implemented
	if (NumTagSets() > 0) {
		cout << "\n NodeManagerT::GetUnknowns: not implemented for XDOF tags" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* check */
	int num_eq = NumEquations(group);
	if (unknowns.Length() != num_eq) throw ExceptionT::kGeneralFail;

	/* loop over groups */
	int checksum = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
		{
			/* field data */
			FieldT& field = *(fFields[i]);
			
			/* check order of field */
			if (field.Order() < order) {
				cout << "\n NodeManagerT::GetUnknowns: order " << order 
				     << " is out of range {0," << field.Order() << "}" << endl;
				throw ExceptionT::kOutOfRange;
			}
			
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
	if (checksum != num_eq) throw ExceptionT::kGeneralFail;
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
	if (source.Length() != target.Length()) ExceptionT::SizeMismatch();

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
	if (fCurrentCoords)
		fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);	
}

#if 0
/* duplicate nodes */
void NodeManagerT::DuplicateNodes(const iArrayT& nodes, 
	iArrayT& new_node_tags)
{
//TEMP
cout << "\n NodeManagerT::DuplicateNodes: not implemented" << endl;
throw ExceptionT::kGeneralFail;
#pragma unused(nodes)
#pragma unused(new_node_tags)

	/* check */
	if (nodes.Length() != new_node_tags.Length()) throw ExceptionT::kSizeMismatch;
	int InitialLength;
	
	/* generate the nodes */
	ModelManagerT* model = fFEManager.ModelManager ();
	model->DuplicateNodes (nodes, new_node_tags, fNumNodes);

	/* set new array lengths */
	fNumNodes += nodes.Length() ;
	fEqnos.Resize(fNumNodes);	
	fDisp.Resize(fNumNodes);

	/* copy data from old nodes to duplicated nodes */
	for( int i = 0 ; i < nodes.Length() ; i++ ) 
	{
		fEqnos.SetRow(new_node_tags[i],0) ;		
		fDisp.CopyRowFromRow(new_node_tags[i],nodes[i]) ;
	}
	
	/* check if initial conditions need duplication */
	InitialLength = fIC.Length() ;
	for( int i = 0 ; i < InitialLength ; i++ ) 
		for( int j = 0 ; j < nodes.Length() ; j++ ) 
		{
			if( fIC[i].Node() == nodes[j] )
			{
				fIC.Resize(fIC.Length()+1) ;
				fIC.Last() = fIC[i] ;
			}
		}
	
	/* check if kinematic boundary conditions need duplication */
	InitialLength = fKBC.Length() ;
	for( int i = 0 ; i < InitialLength ; i++ ) 
		for( int j = 0 ; j < nodes.Length() ; j++ ) 
		{
			if( fKBC[i].Node() == nodes[j] )
			{
				fKBC.Resize(fKBC.Length()+1) ;
				fKBC.Last() = fKBC[i] ;
			}
		}

	/* check if force boundary conditions need duplication */
	InitialLength = fFBC.Length() ;
	for( int i = 0 ; i < InitialLength ; i++ ) 
		for( int j = 0 ; j < nodes.Length() ; j++ ) 
		{
			if( fFBC[i].Node() == nodes[j] )
			{
				fFBC[i].SplitForce();
				fFBC.Resize(fFBC.Length()+1) ;
				fFBC.Last() = fFBC[i] ;
				fFBCValues.Resize(fFBC.Length()); 
				fFBCEqnos.Resize(fFBC.Length());
				fFBCEqnosSet = 0;
			}
		}

	/* (re-)set start tag for external DOF */
	XDOF_ManagerT::SetStartTag(fNumNodes);

	/* resize averaging work space */
	SetNumAverageRows(fNumNodes);
}
#endif

//TEMP - trap parallel execution with XDOF
void NodeManagerT::XDOF_Register(DOFElementT* group, const iArrayT& numDOF)
{
	//TEMP - parallel execution not yet supported
	if (fFEManager.Size() > 1)
	{
		cout << "\n NodeManagerT::NodeManagerT: not for parallel execution" << endl;
		throw ExceptionT::kGeneralFail;
	}
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
	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0) {
		cout << "\n NodeManagerT::XDOF_SetLocalEqnos: group has not fields: " 
		     << group << endl;
		throw ExceptionT::kGeneralFail;
	}

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
				{
					cout << "\n NodeManagerT::XDOF_SetLocalEqnos: could not resolve tag into set: " 
					     << tag << endl;
					throw ExceptionT::kGeneralFail;
				}

				/* equations from tag set */
				eqnos_source = fXDOF_Eqnos[tag_set];
			}		
			
			/* dimension */
			int ndof = eqnos_source->MinorDim();
			eq_count += ndof;

			/* check number of assigned equations */
			if (eq_count > neq)
			{
				cout << "\n NodeManagerT::XDOF_SetLocalEqnos: error assigning equations" << endl;
				throw ExceptionT::kSizeMismatch;
			}
		
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
	/* check */
	if (nodes.MajorDim() != eqnos.MajorDim()) throw ExceptionT::kSizeMismatch;

	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0) {
		cout << "\n NodeManagerT::XDOF_SetLocalEqnos: group has not fields: " 
		     << group << endl;
		throw ExceptionT::kGeneralFail;
	}

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
					{
						cout << "\n NodeManagerT::XDOF_SetLocalEqnos: could not resolve tag into set: " 
						     << tag << endl;
						throw ExceptionT::kGeneralFail;
					}
				
					/* equations from tag set */
					eqnos_source = fXDOF_Eqnos[tag_set];
				}

				/* dimension */
				int ndof = eqnos_source->MinorDim();
				eq_count += ndof;

				/* check number of assigned equations */
				if (eq_count > neq)
				{
					cout << "\n NodeManagerT::XDOF_SetLocalEqnos: error assigning equations" << endl;
					throw ExceptionT::kSizeMismatch;
				}

				/* copy equations */
				eqnos_source->RowCopy(tag - tag_offset, peq);
			
				/* next */
				peq += ndof;
			}
		}
		
		/* check */
		if (eq_count != neq) {
			cout << "\n NodeManagerT::XDOF_SetLocalEqnos: the number of assigned equations " << eq_count 
			     << "\n     is not equal to the number of element equations " << neq << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

/* collection equation numbers for mixed connectivities. */
void NodeManagerT::XDOF_SetLocalEqnos(int group, const RaggedArray2DT<int>& nodes, 
	RaggedArray2DT<int>& eqnos) const
{
	/* check */
	if (nodes.MajorDim() != eqnos.MajorDim()) throw ExceptionT::kSizeMismatch;

	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0) {
		cout << "\n NodeManagerT::XDOF_SetLocalEqnos: group has not fields: " 
		     << group << endl;
		throw ExceptionT::kGeneralFail;
	}

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
					{
						cout << "\n NodeManagerT::XDOF_SetLocalEqnos: could not resolve tag into set: " 
						     << tag << endl;
						throw ExceptionT::kGeneralFail;
					}
				
					/* equations from tag set */
					eqnos_source = fXDOF_Eqnos[tag_set];
				}

				/* dimension */
				int ndof = eqnos_source->MinorDim();
				eq_count += ndof;

				/* check number of assigned equations */
				if (eq_count > neq)
				{	
					cout << "\n NodeManagerT::XDOF_SetLocalEqnos: error assigning equations" << endl;
					throw ExceptionT::kSizeMismatch;
				}

				/* copy equations */
				eqnos_source->RowCopy(tag - tag_offset, peq);
				
				/* next */
				peq += ndof;
			}
		}

		/* check */
		if (eq_count != neq) {
			cout << "\n NodeManagerT::XDOF_SetLocalEqnos: the number of assigned equations " << eq_count 
			     << "\n     is not equal to the number of element equations " << neq << endl;
			throw ExceptionT::kGeneralFail;
		}		
	}
}

/* information about subordinate parameter lists */
void NodeManagerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the fields */
	sub_list.AddSub("field", ParameterListT::OnePlus);
	
	/* list of history node ID's */
	sub_list.AddSub("history_node_ID", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NodeManagerT::NewSub(const StringT& list_name) const
{
	if (list_name == "field")
		return new FieldT(fFieldSupport);
	else if (list_name == "history_node_ID")
		return new IntegerListT("history_node_ID");
	else
		return ParameterInterfaceT::NewSub(list_name);
}

/**********************************************************************
 * Protected
 **********************************************************************/

void NodeManagerT::EchoCoordinates(ifstreamT& in, ostream& out)
{
	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* read coordinates */
	if (model->DatabaseFormat() == IOBaseT::kTahoe)
		model->ReadInlineCoordinates (in);
	else
		model->ReadCoordinates();
		
	/* set pointer */	
	fInitCoords = &(model->Coordinates());

	/* check element groups to see if node data should be 
	   adjusted to be 2D, some element groups require
	   fNumSD == fNumDOF */
	if (NumSD() == 3 && model->AreElements2D())
	{
		out << " Number of nodal points. . . . . . . . . . . . . = " << NumNodes() << '\n';
		out << " Number of spatial dimensions. . . . . . . . . . = " << NumSD() << '\n';
		out << "\n Adjusting nodal data to 2D" << endl;
		cout << "\n NodeManagerT::EchoCoordinates: WARNING: Adjusting nodal data to 2D" << endl;
		model->AdjustCoordinatesto2D();
	}
	  
	/* print main header */
	out << " Number of nodal points. . . . . . . . . . . . . = " << NumNodes() << '\n';
	out << " Number of spatial dimensions. . . . . . . . . . = " << NumSD() << '\n';

	/* verbose output */
	if (fFEManager.PrintInput())
	{
		int d_width = out.precision() + kDoubleExtra;
	
		/* write header */
		out << setw(kIntWidth) << "node"
		    << setw(kIntWidth) << "gl.node"
		    << setw(kIntWidth) << "proc";		
		for (int i = 0; i < NumSD(); i++)
			out << setw(d_width - 2) << "x[" << i + 1 << "]";
		out << '\n';

		/* arrays */
		const ArrayT<int>* processor = ProcessorMap();
		const dArray2DT& init_coords = InitialCoordinates();
		const ArrayT<int>* node_map = fFEManager.NodeMap();
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

/* echo field data */
void NodeManagerT::EchoFields(ifstreamT& in, ostream& out)
{
	const char caller[] = "NodeManagerT::EchoFields";

	/* no predefined fields */
	if (fFEManager.Analysis() == GlobalT::kMultiField)
	{
		/* number of fields */
		int num_fields = -1;
		in >> num_fields;

		/* initialize list */
		fFields.Dimension(num_fields);
		fFields = NULL;
		fMessageID.Dimension(num_fields);
		
		/* loop over fields */
		for (int i = 0; i < fFields.Length(); i++)
		{
			int index = -1, group_num = -1, cont_num = -1, ndof = -1;
			StringT name;
			
			/* read: field index */
			in >> index;
			
			/* read: [name] [group] [controller] [ndof] */
			in >> name >> group_num >> cont_num >> ndof;
			index--;
			group_num--;
			cont_num--;
			if (fFields[index] != NULL)
				ExceptionT::BadInputValue(caller, "field at index %d is already set", index+1);

			/* check for field with same name */
			if (Field(name)) 
				ExceptionT::BadInputValue(caller, "field name %s already exists", name.Pointer());

			/* read: dof labels */
			ArrayT<StringT> labels(ndof);
			for (int j = 0; j < labels.Length(); j++)
				in >> labels[j];
			
			/* get integrator */
			const nIntegratorT* controller = fFEManager.nIntegrator(cont_num);
			if (!controller) ExceptionT::GeneralFail(caller);

			/* new field */			
			FieldT* field = new FieldT(fFieldSupport);
			field->Initialize(name, ndof, *controller);
			field->SetLabels(labels);
			field->SetGroup(group_num);
			field->Dimension(NumNodes(), false);
			field->Clear();
			field->WriteParameters(out);

			/* coordinate update field */
			if (name == "displacement") {
				if (fCoordUpdate) ExceptionT::BadInputValue(caller, "\"displacement\" field already set");
				fCoordUpdate = field;
				fCurrentCoords = new dArray2DT;
				fCurrentCoords_man.SetWard(0, *fCurrentCoords, NumSD());
				fCurrentCoords_man.SetMajorDimension(NumNodes(), false);
				(*fCurrentCoords) = InitialCoordinates();
			}

			/* echo initial/boundary conditions */
			EchoInitialConditions(*field, in, out);
			EchoKinematicBC(*field, in, out);
			EchoKinematicBCControllers(*field, in, out);
			EchoForceBC(*field, in, out);
			EchoForceBCControllers(*field, in, out);

			/* store */
			fFields[index] = field;
			
			/* set up communication of field */
			fMessageID[index] = fCommManager.Init_AllGather(fFields[index]->Update());
		}
	}
	else /* legacy code - a single predefined integrator */
	{
		/* just one */
		fFields.Dimension(1);
		fFields = NULL;
		fMessageID.Dimension(1);
		const nIntegratorT* controller = fFEManager.nIntegrator(0);
		if (!controller) ExceptionT::GeneralFail(caller);
		
		/* field config set by analysis type */
		FieldT* field = NULL;
		switch(fFEManager.Analysis())
		{
			case GlobalT::kLinDynamic:
			case GlobalT::kNLDynamic:
			case GlobalT::kLinExpDynamic:
			case GlobalT::kNLExpDynamic:
			case GlobalT::kNLExpDynKfield:
			case GlobalT::kNLStatic:
			case GlobalT::kNLStaticKfield:
			case GlobalT::kLinStatic:
			{
				field = new FieldT(fFieldSupport);
				field->Initialize("displacement", NumSD(), *controller);

				/* label list */
				ArrayT<StringT> labels(field->NumDOF());
				const char* disp[] = {"D_X", "D_Y", "D_Z"};
				for (int i = 0; i < labels.Length(); i++)
					labels[i] = disp[i];

				/* set labels */
				field->SetLabels(labels);				
				break;
			}
			case GlobalT::kLinTransHeat:
			case GlobalT::kLinStaticHeat:
			case GlobalT::kNLStaticHeat:
			case GlobalT::kNLTransHeat:
			{
				field = new FieldT(fFieldSupport);
				field->Initialize("temperature", 1, *controller);

				/* set labels */
				ArrayT<StringT> labels(1);
				labels[0] = "T";
				field->SetLabels(labels);
				break;
			}
			default:
				ExceptionT::GeneralFail(caller, "unrecognized analysis code %d", fFEManager.Analysis());
		}
		
		/* check */
		if (!field) throw ExceptionT::kGeneralFail;
		
		/* just one group */
		field->SetGroup(0);
		
		/* dimension */
		field->Dimension(NumNodes(), false);
		field->Clear();
		
		/* clear all equation numbers */
		field->Equations() = FieldT::kInit;

		/* coordinate update field */
		if (field->Name() == "displacement") {
			if (fCoordUpdate) ExceptionT::BadInputValue(caller, "\"displacement\" field already set");
			fCoordUpdate = field;
			fCurrentCoords = new dArray2DT(InitialCoordinates());
		}

		/* initial/boundary conditions */
		EchoInitialConditions(*field, in, out);
		EchoKinematicBC(*field, in, out);
		EchoKinematicBCControllers(*field, in, out);
		EchoForceBC(*field, in, out);
		EchoForceBCControllers(*field, in, out);
		
		/* store */
		fFields[0] = field;

		/* set up communication of field */
		fMessageID[0] = fCommManager.Init_AllGather(fFields[0]->Update());
	}
}

void NodeManagerT::EchoInitialConditions(FieldT& field, ifstreamT& in, ostream& out)
{
	out << "\n Initial conditions:\n\n";

	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* initial condition cards */
	ArrayT<IC_CardT>& cards = field.InitialConditions();

	/* account for text file name instead of data */
	ifstreamT tmp;
	ifstreamT& in2 = model->OpenExternal(in, tmp, out, true, "NodeManagerT::EchoInitialConditions: could not open file");
	
	/* read */
	int numIC_sets;
	in2 >> numIC_sets;

	if (numIC_sets > 0)
	{
		/* allocate temp space */
		ArrayT<iArrayT> nodes (numIC_sets);
		iArray2DT data (numIC_sets, 2);
		dArrayT values (numIC_sets);

		int count = model->ReadCards(in2, out, nodes, data, values);
		cards.Dimension(count);

	    /* set cards */
	    for (int i=0, dex=0; i < numIC_sets; i++)
		if (nodes[i].Length() > 0) {
			int *pn = nodes[i].Pointer();
			
			/* ID == -2 (-1 shifted by 1) means "all" */
			if (*pn == -2) //offset
			{
				/* resize array */
				nodes[i].Dimension(NumNodes());
				nodes[i].SetValueToPosition();
				pn = nodes[i].Pointer();
				
				/* need more cards */
				cards.Resize(cards.Length() - 1 + nodes[i].Length());
			}
			
			int dof = data (i, 0) - 1; // offset
			int order = data (i, 1) - 1; // KBC code -> derivative order
			for (int j=0; j < nodes[i].Length(); j++)
				cards[dex++].SetValues(*pn++, dof, order, values[i]);
			}
	  }

	/* echo */
	out << " Number of initial condition cards . . . . . . . = ";
	out << cards.Length() << "\n\n";
	if (fFEManager.PrintInput() && cards.Length() > 0)
	{
		IC_CardT::WriteHeader(out);
		for (int i = 0; i < cards.Length(); i++)
			cards[i].WriteValues(out);
	}
}

void NodeManagerT::EchoKinematicBC(FieldT& field, ifstreamT& in, ostream& out)
{
	out << "\n Kinematic boundary conditions:\n\n";

	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* kinematic boundary condition cards */
	ArrayT<KBC_CardT>& cards = field.KinematicBC();

	/* account for text file name instead of data */
	ifstreamT tmp;
	ifstreamT& in2 = model->OpenExternal(in, tmp, out, true, "NodeManagerT::EchoKinematicBC: could not open file");

	int numKBCsets;
	in2 >> numKBCsets;
	if (numKBCsets > 0)
	  {
	    /* allocate temp space */
	    ArrayT<iArrayT> nodes (numKBCsets);
	    iArray2DT data (numKBCsets, 3);
	    dArrayT values (numKBCsets);

	    int count = model->ReadCards (in2, out, nodes, data, values);
	    cards.Dimension(count);

	    /* set cards */
	    for (int i=0, dex=0; i < numKBCsets; i++)
		if (nodes[i].Length() > 0) {
			int *pn = nodes[i].Pointer();
			
			/* ID == -2 (-1 shifted by 1) means "all" */
			if (*pn == -2) //offset
			{
				/* resize array */
				nodes[i].Dimension(NumNodes());
				nodes[i].SetValueToPosition();
				pn = nodes[i].Pointer();

				/* need more cards - only counted an ID of -1 */
				cards.Resize(cards.Length() - 1 + nodes[i].Length());
			}

			int dof = data (i, 0) - 1; // offset
			KBC_CardT::CodeT code = cards[dex].int_to_CodeT (data (i, 1));
			int ltf = data (i, 2) - 1; // offset
			for (int j=0; j < nodes[i].Length(); j++)
				cards[dex++].SetValues (*pn++, dof, code, ltf, values[i]);
		}
	  }

	/* finalize BC's/echo data */
	out << "\n Number of kinematic boundary condition cards. . = ";
	out << cards.Length() << "\n\n";
	if (fFEManager.PrintInput() && cards.Length() > 0) 
		KBC_CardT::WriteHeader(out);

	for (int i = 0; i < cards.Length(); i++)
	{
		/* card */
		KBC_CardT& card = cards[i];

		/* set schedule pointer */
		const ScheduleT* schedule = (card.Code() != KBC_CardT::kFix) ? Schedule(card.ScheduleNum()) : NULL;
		card.SetSchedule(schedule);
		
		/* echo values */
		if (fFEManager.PrintInput()) card.WriteValues(out);
	}

	/* check node numbers */
	for (int i = 0; i < cards.Length(); i++)
		if (cards[i].Node() < 0 || cards[i].Node() >= NumNodes()) {
			cout << "\n NodeManagerT::EchoKinematicBC: node number is out of range: " 
			     << cards[i].Node() + 1 << endl;
			throw ExceptionT::kOutOfRange;
		}
}

void NodeManagerT::EchoForceBC(FieldT& field, ifstreamT& in, ostream& out)
{
	const char caller[] = "NodeManagerT::EchoForceBC";

	/* print header */
	out << "\n Force boundary conditions:\n\n";

	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* nodal force cards */
	ArrayT<FBC_CardT>& cards = field.ForceBC();

	/* account for text file name instead of data */
	ifstreamT tmp;
	ifstreamT& in2 = model->OpenExternal(in, tmp, out, true, "NodeManagerT::EchoKinematicBC: could not open file");

	int numFBCsets;
	in2 >> numFBCsets;

	/* check */
	if (numFBCsets < 0)
		ExceptionT::BadInputValue(caller, "expecting numFBC_sets > 0: %d", numFBCsets);

	if (numFBCsets > 0)
	  {
	    /* allocate temp space */
	    ArrayT<iArrayT> nodes (numFBCsets);
	    iArray2DT data (numFBCsets, 2);
	    dArrayT values (numFBCsets);

	    int count = model->ReadCards(in2, out, nodes, data, values);
		cards.Dimension (count);

	    /* set cards */
	    for (int i=0, dex=0; i < numFBCsets; i++)
		if (nodes[i].Length() > 0) {
			int *pn = nodes[i].Pointer();
		
			/* ID == -2 (-1 shifted by 1) means "all" */
			if (*pn == -2) //offset
			{
				/* resize array */
				nodes[i].Dimension(NumNodes());
				nodes[i].SetValueToPosition();
				pn = nodes[i].Pointer();

				/* need more cards - only counted an ID of -1 */
				cards.Resize(cards.Length() - 1 + nodes[i].Length());
			}

			int dof = data (i, 0) - 1; //offset
			int ltf = data (i, 1) - 1; //offset

			for (int j=0; j < nodes[i].Length(); j++)
				cards[dex++].SetValues (*this, *pn++, dof, ltf, values[i]);
		}
	  }

	/* echo */
	out << " Number of nodal force cards . . . . . . . . . . = ";
	out << cards.Length() << "\n\n";
	if (fFEManager.PrintInput() && cards.Length() > 0)
	{
		FBC_CardT::WriteHeader(out);
		for (int i = 0; i < cards.Length(); i++)
			cards[i].WriteValues(out);
	}

	/* check node numbers */
	for (int i = 0; i < cards.Length(); i++)
		if (cards[i].Node() < 0 || cards[i].Node() >= NumNodes()) 
			ExceptionT::OutOfRange(caller, "node number is out of range: %d", cards[i].Node() + 1);
}

void NodeManagerT::EchoHistoryNodes(ifstreamT& in, ostream &out)
{
	out << "\n N o d a l   H i s t o r i e s :\n\n";

	/* model database */
	ModelManagerT* model = fFEManager.ModelManager();
	
	/* list of external nodes */
	const ArrayT<int>* p_ex_nodes = fCommManager.ExternalNodes();
	iArrayT ex_nodes;
	if (p_ex_nodes) ex_nodes.Alias(*p_ex_nodes);

	/* read node set indexes */
	ArrayT<StringT> node_set_ID;
	model->NodeSetList(in, node_set_ID);

	/* just read node sets. output of history information is
	 * configured in RegisterOutput */
	int num_sets = node_set_ID.Length();
	fHistoryNodeSetIDs.Dimension(num_sets);
	for (int i = 0; i < num_sets; i++)
	{
		/* read node set */
		const iArrayT& node_set = model->NodeSet(node_set_ID[i]);
		
		/* slow-but-steady way */
		ArrayT<bool> is_external(node_set.Length());
		is_external = true;
		int count = 0;
		for (int j = 0; j < node_set.Length(); j++)
			if (!ex_nodes.HasValue(node_set[j]))
			{
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
			ExceptionT::BadInputValue("NodeManagerT::EchoHistoryNodes", 
				"error initializing node set %d as model element ID %d", node_set_ID[i].Pointer(), ID.Pointer());

		fHistoryNodeSetIDs[i] = ID;
	}
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
	switch(code)
	{
		case KBC_ControllerT::kK_Field:
			return new K_FieldT(*this);

		case KBC_ControllerT::kBimaterialK_Field:	
			return new BimaterialK_FieldT(*this);

		case KBC_ControllerT::kMappedPeriodic:	
			return new MappedPeriodicT(*this, field);

		case KBC_ControllerT::kTiedNodes:
		{
			TiedNodesT* kbc = new TiedNodesT(*this, field);
			return kbc;
		}
		case KBC_ControllerT::kPeriodicNodes:
		{
			PeriodicNodesT* kbc = new PeriodicNodesT(*this, field);
			return kbc;
		}
		case KBC_ControllerT::kScaledVelocityNodes:
		{
			ScaledVelocityNodesT* kbc = new ScaledVelocityNodesT(*this, field);
			return kbc;
		}
		case KBC_ControllerT::kSetOfNodesKBC:
		{
			SetOfNodesKBCT* kbc = new SetOfNodesKBCT(*this, field);
			return kbc;
		}
		case KBC_ControllerT::kTorsion:
		{
			TorsionKBCT* kbc = new TorsionKBCT(*this, fFEManager.Time());
			return kbc;
		}
		case KBC_ControllerT::kConyevor:
		{
			ConveyorT* kbc = new ConveyorT(*this, field);
			return kbc;
		}
		default:
			ExceptionT::BadInputValue("NodeManagerT::NewKBC_Controller", 
				"KBC controller code %d is not supported", code);
	}
	return NULL;
}

FBC_ControllerT* NodeManagerT::NewFBC_Controller(FieldT& field, int code)
{
  const char caller[] = "NodeManagerT::NewFBC_Controller";

	/* displacement field */
	const dArray2DT& disp = field[0];

	/* velocity field */
	dArray2DT* velocity = NULL;
	if (field.Order() > 0)
		velocity = &(field[1]);
		
	/* const equations array */
	const iArray2DT& eqnos = field.Equations();
	
	/* current coordinates */
	const dArray2DT& coords = CurrentCoordinates();

	FBC_ControllerT* fbc = NULL;
	switch(code)
	{
		case FBC_ControllerT::kPenaltyWall:
			fbc = new PenaltyWallT(fFEManager, field.Group(), eqnos, coords, disp, velocity);
			break;

		case FBC_ControllerT::kAugLagWall:
			fbc = new AugLagWallT(fFEManager, this, field, coords, disp);
			break;

		case FBC_ControllerT::kPenaltySphere:	
			fbc = new PenaltySphereT(fFEManager, field.Group(), eqnos, coords, disp, velocity);
			break;

		case FBC_ControllerT::kPenaltyCylinder:	
			fbc = new PenaltyCylinderT(fFEManager, field.Group(), eqnos, coords, disp, velocity);
			break;

		case FBC_ControllerT::kAugLagSphere:	
			fbc = new AugLagSphereT(fFEManager, this, field, coords, disp);
			break;

		case FBC_ControllerT::kMFPenaltySphere:	
			fbc = new MFPenaltySphereT(fFEManager, field.Group(), eqnos, coords, disp, velocity);
			break;
	    case FBC_ControllerT::kMFAugLagMult:
	    	fbc = new MFAugLagMultT(fFEManager, this, field, coords, disp);
	    	break;
		default:
			ExceptionT::BadInputValue(caller, "FBC controller code %d is not supported", code);
	}
	
	/* set time integrator */
	if (fbc) {
		const eIntegratorT& e_integrator = field.nIntegrator().eIntegrator();
		fbc->SetController(&e_integrator);
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

void NodeManagerT::EchoKinematicBCControllers(FieldT& field, ifstreamT& in, ostream& out)
{	
	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* account for text file name instead of data */
	ifstreamT tmp;
	ifstreamT& in2 = model->OpenExternal(in, tmp, out, true, "NodeManagerT::EchoKinematicBCControllers: could not open file");

	int numKBC;
	in2 >> numKBC;

	/* print header */
	out << "\n K i n e m a t i c   B C   C o n t r o l l e r s :\n\n";
	out << " Number of controllers . . . . . . . . . . . . . = ";
	out << numKBC << '\n' << endl;

	/* construct */
	for (int i = 0; i < numKBC; i++)
	{
		int num, KBC_type;
		in2 >> num >> KBC_type; num--;
		
		/* construct */
		KBC_ControllerT* controller = NewKBC_Controller(field, KBC_type);
		
		/* initialize */
		controller->Initialize(in2);
		
		/* store */
		field.AddKBCController(controller);
		
		/* special handling of conveyor */
		if (false && KBC_type == KBC_ControllerT::kConyevor) {
		
			/* force all geometry to be changing */
			fFEManager.OutputManager()->SetChangingFlag(IOManager::kForceChanging);
		}
	}

	/* echo parameters */
	const ArrayT<KBC_ControllerT*>& controllers = field.KBC_Controllers();
	for (int j = 0; j < numKBC; j++)
	{
		out << " Controller: " << j+1 << '\n';
	
		/* write parameters */
		controllers[j]->WriteParameters(out);
		out << '\n';
	}

	/* flush output */
	out.flush();
}

void NodeManagerT::EchoForceBCControllers(FieldT& field, ifstreamT& in, ostream& out)
{
	/* model manager */
	ModelManagerT* model = fFEManager.ModelManager();

	/* account for text file name instead of data */
	ifstreamT tmp;
	ifstreamT& in2 = model->OpenExternal(in, tmp, out, true, "NodeManagerT::EchoForceBCControllers: could not open file");
	
	/* FBC controllers */
	int numFBCcont;
	in2 >> numFBCcont;

	/* print header */
	out << "\n N o d a l   F o r c e   C o n t r o l l e r s :\n\n";
	out << " Number of controllers . . . . . . . . . . . . . = ";
	out << numFBCcont << '\n' << endl;

	for (int i = 0; i < numFBCcont; i++)
	{
		int num, type;
		in2 >> num >> type; num--;
		
		/* construct */
		FBC_ControllerT* controller = NewFBC_Controller(field, type);
		
		/* echo data */
		controller->EchoData(in2, out);

		/* initialize */
		controller->Initialize();
		
		/* store */
		field.AddFBCController(controller);
	}

	/* flush output */
	out.flush();
}

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
