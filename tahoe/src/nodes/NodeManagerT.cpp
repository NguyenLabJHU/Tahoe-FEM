/* $Id: NodeManagerT.cpp,v 1.19 2003-01-21 16:54:17 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "NodeManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <limits.h>
#include <ctype.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "nControllerT.h"
#include "eControllerT.h"
#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "PartitionT.h"
#include "ReLabellerT.h"
#include "OutputSetT.h"

#include "FieldSupportT.h"
#include "FieldT.h"

/* force BC controllers */
#include "AugLagWallT.h"
#include "PenaltyWallT.h"
#include "PenaltySphereT.h"
#include "AugLagSphereT.h"
#include "MFPenaltySphereT.h"

/* kinematic BC controllers */
#include "K_FieldT.h"
#include "BimaterialK_FieldT.h"
#include "MappedPeriodicT.h"
#include "TiedNodesT.h"
#include "SymmetricNodesT.h"
#include "PeriodicNodesT.h"

using namespace Tahoe;

/* constructor */
NodeManagerT::NodeManagerT(FEManagerT& fe_manager):
	fFEManager(fe_manager),
	fInitCoords(NULL),
	fCoordUpdate(NULL),
	fCurrentCoords(NULL)
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
	int neq = 0;
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			neq += fFields[i]->NumActiveEquations();
			
	/* get XDOF equations */
	neq += XDOF_ManagerT::NumEquations(group);
				
	return neq; 
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

	/* external nodes (parallel execution) */
	EchoExternalNodes(out);

	/* history nodes */
	EchoHistoryNodes(in, out);
}

/* register data for output */
void NodeManagerT::RegisterOutput(void)
{
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
			/* field labels */
			const ArrayT<StringT>& labels = fFields[j]->Labels();
			int ndof = labels.Length();
			
			/* output labels */
			ArrayT<StringT> n_labels(2*ndof); /* 2x : the DOF's and their conjugates */
			for (int i = 0; i < ndof; i++)
			{
				n_labels[i] = labels[i]; /* DOF */
				n_labels[i+ndof].Append("F_", n_labels[i]); /* conjugate */
			}

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
					model.UpdateConnectivity(ID, new_conn, false);
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
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* from fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->Connectivities(connects_1, connects_2);
}

/* return the implicit-explicit flag for the given group */
ControllerT::ImpExpFlagT NodeManagerT::ImplicitExplicit(int group) const
{
	ControllerT::ImpExpFlagT flag = ControllerT::kExplicit;

	/* loop over fields in the group */
	for (int i = 0; flag == ControllerT::kExplicit && i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			flag = fFields[i]->nController().ImplicitExplicit();

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
			break;					
		}
		default:
			cout << "\n NodeManagerT::RegisterCoordinates: not a coordinate type: " 
			     << array.Type() << endl;
			throw ExceptionT::kGeneralFail;
	}
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
	{
		/* should be allocated */
		if (!fCurrentCoords)
			ExceptionT::GeneralFail("NodeManagerT::InitStep", "current coords not initialized");
	
		/* update */
		fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
	}	
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
		FieldSupportT support(fFEManager);
		for (int i = 0; i < fFields.Length(); i++)
			if (fFields[i]->Group() == group)
				fFields[i]->FormLHS(support, sys_type);
	}
}
	
/* compute the nodal contribution to the residual force vector */
void NodeManagerT::FormRHS(int group)
{
	FieldSupportT support(fFEManager);
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->FormRHS(support);
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
	/* active equation numbers */
	int eq_start = fFEManager.ActiveEquationStart(group);
	int num_eq   = NumEquations(group);
	
	/* update fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->Update(update, eq_start, num_eq);

	/* external nodes */
	if (Size() > 1)
		for (int i = 0; i < fFields.Length(); i++)
		{
			FieldT& field = *(fFields[i]);
		
			/* field data */
			nControllerT& integrator = field.nController();
			dArray2DT& update = field.ExternalUpdate();
			iArray2DT& xeqnos = field.ExternalEquations();
		
			/* send update for outgoing */
			fFEManager.SendExternalData(integrator.MappedCorrectorField(field));
	
			/* receive update for incoming */
			fFEManager.RecvExternalData(update);

			/* apply update - external nodes to be updated are marked with 1 */
			integrator.MappedCorrector(field, fExNodes, xeqnos, update);
		}

	/* update current configurations */
	if (fCoordUpdate && fCoordUpdate->Group() == group)
	{
		/* should be allocated */
		if (!fCurrentCoords) {
			cout << "\n NodeManagerT::Update: current coords not initialized" << endl;
			throw ExceptionT::kGeneralFail;
		}
	
		/* update */
		fCurrentCoords->SumOf(InitialCoordinates(), (*fCoordUpdate)[0]);
	}	
	
	/* inherited - update external DOF */
	XDOF_ManagerT::Update(group, update);
}

/* update history */
void NodeManagerT::CloseStep(int group)
{
	/* output nodal history */
//NOTE - moved to WriteOutput
//	WriteNodalHistory(group);

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
		fFields[i]->InitialCondition();
}

void NodeManagerT::ReadRestart(ifstreamT& in)
{
	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->ReadRestart(in);
}

void NodeManagerT::WriteRestart(ofstreamT& out) const
{
	/* apply to fields */
	for (int i = 0; i < fFields.Length(); i++)
		fFields[i]->WriteRestart(out);
}

/* reset displacements (and configuration to the last known solution) */
void NodeManagerT::ResetStep(int group)
{
	/* reset fields */
	for (int i = 0; i < fFields.Length(); i++)
		if (fFields[i]->Group() == group)
			fFields[i]->ResetStep();

	/* inherited - reset external DOF */
	XDOF_ManagerT::Reset(group);
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

	/* nodal histories */
//NOTE - moved here from CloseStep	
	WriteNodalHistory();
}

void NodeManagerT::SetEquationNumbers(int group)
{
	/* collect fields in the group */
	ArrayT<FieldT*> fields;
	CollectFields(group, fields);
	if (fields.Length() == 0) {
		cout << "\n NodeManagerT::SetEquationNumbers: group has not fields: " 
		     << group << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* initialize equations numbers arrays */
	ArrayT<iArray2DT> group_eqnos(fields.Length());
	for (int i = 0; i < fields.Length(); i++)
		fields[i]->InitEquations(group_eqnos[i]);

	/* mark external nodes */
	if (fExNodes.Length() > 0)
		for (int i = 0; i < fields.Length(); i++)
		{
			/* field data */
			iArray2DT&  eqnos = group_eqnos[i];
			iArray2DT& xeqnos = fields[i]->ExternalEquations();
		
			/* collect equation numbers */
			xeqnos.RowCollect(fExNodes, eqnos);
	
			/* mark active external nodes for update of active
			 * external using controller */
			for (int j = 0; j < xeqnos.Length(); j++)
				if (xeqnos[j] >= FieldT::kInit) 
					xeqnos[j] = 1;

			/* mark all external as inactive for setting local
			 * equation numbers */
			for (int j = 0; j < fExNodes.Length(); j++)
				eqnos.SetRow(fExNodes[j], FieldT::kExternal);	
		}

	/* assign active equation numbers node-by-node across fields
	 * in the group */
	iArrayT num_active(fields.Length());
	num_active = 0;
	int num_eq = 0;
	int nnd = NumNodes();
	for (int i = 0; i < nnd; i++)
		for (int j = 0; j < fields.Length(); j++)
		{
			int& count = num_active[j];
			int  ndof = fields[j]->NumDOF();
			int* peq = (group_eqnos[j])(i);
			for (int k = 0; k < ndof; k++)
			{
				/* active equation */
				if (*peq >= FieldT::kInit) 
				{
					*peq = ++num_eq;
					count++;
				}
				peq++;
			}
		}

	/* set equation arrays */
	for (int i = 0; i < fields.Length(); i++)
		fields[i]->SetEquations(group_eqnos[i], num_active[i]);
	
	/* assign equation numbers to XDOF's */
	XDOF_ManagerT::SetEquations(group, num_eq);
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
	{
		cout << "\n NodeManagerT::RenumberEquations: expecting to renumber "
		     << NumEquations(group) << '\n'
		     <<   "     equations, but hit " << numtest << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* rearrange equations if needed */
	CheckEquationNumbers(group);

	cout << "\n NodeManagerT::RenumberEquations: done" << endl;
}

void NodeManagerT::SetEquationNumberScope(int group, GlobalT::EquationNumberScopeT scope)
{
	//TEMP - external DOF's no tested with other scopes
	if (scope != GlobalT::kLocal && fExNodes.Length() > 0 && NumTagSets() > 0)
	{
		cout << "\n NodeManagerT::SetEquationNumberScope: external DOF only\n"
		     <<   "     verified with local numbering" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* switch numbering scope - with external nodes */
	if (scope == GlobalT::kGlobal && fExNodes.Length() > 0)
	{
		/* shift local equation numbers */
		int start = fFEManager.GlobalEquationStart(group);
		int shift = start - 1;
		
		/* loop over fields */
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

				/* collect external equation numbers */
				iArray2DT external_equations(fExNodes.Length(), fFields[j]->NumDOF());
				fFEManager.SendRecvExternalData(eqnos, external_equations);

				/* write-in */
				eqnos.Assemble(fExNodes, external_equations);
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

	int* ien = nodes.Pointer();
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

	int* ien = nodes.Pointer();
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

	int* ien = nodes.Pointer();
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

	/* set node to processor map */
	iArrayT full_range(NumNodes());
	full_range.SetValueToPosition();
	fFEManager.NodeToProcessorMap(full_range, fProcessor);
	full_range.Free();

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
		const dArray2DT& init_coords = InitialCoordinates();
		const iArrayT* node_map = fFEManager.NodeMap();
		for (int i = 0; i < init_coords.MajorDim(); i++)
		{
			out << setw(kIntWidth) << i+1
			    << setw(kIntWidth) << ((node_map) ? (*node_map)[i]+1 : i+1)
			    << setw(kIntWidth) << fProcessor[i];		
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
	/* no predefined fields */
	if (fFEManager.Analysis() == GlobalT::kMultiField)
	{
		/* number of fields */
		int num_fields = -1;
		in >> num_fields;

		/* initialize list */
		fFields.Dimension(num_fields);
		fFields = NULL;
		
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
			if (fFields[index] != NULL) {
				cout << "\n NodeManagerT::EchoFields: field at index "
				     << index+1 <<   "     is already set" << endl;
				throw ExceptionT::kBadInputValue; 
			}

			/* check for field with same name */
			if (Field(name)) {
				cout << "\n NodeManagerT::EchoFields: field with name " << name 
				     << " already exists" << endl;
				throw ExceptionT::kBadInputValue;
			}
			
			/* read: dof labels */
			ArrayT<StringT> labels(ndof);
			for (int j = 0; j < labels.Length(); j++)
				in >> labels[j];
			
			/* get integrator */
			nControllerT* controller = fFEManager.nController(cont_num);
			if (!controller) throw ExceptionT::kGeneralFail;

			/* new field */			
			FieldT* field = new FieldT(name, ndof, *controller);
			field->SetLabels(labels);
			field->SetGroup(group_num);
			field->Dimension(NumNodes());
			field->WriteParameters(out);

			/* coordinate update field */
			if (name == "displacement") {
				if (fCoordUpdate) {
					cout << "\n NodeManagerT::EchoFields: \"displacement\" field already set" << endl;
					throw ExceptionT::kBadInputValue;
				}
				fCoordUpdate = field;
				fCurrentCoords = new dArray2DT(InitialCoordinates());
			}

			/* echo initial/boundary conditions */
			EchoInitialConditions(*field, in, out);
			EchoKinematicBC(*field, in, out);
			EchoKinematicBCControllers(*field, in, out);
			EchoForceBC(*field, in, out);
			EchoForceBCControllers(*field, in, out);

			/* store */
			fFields[index] = field;
		}
	}
	else /* legacy code - a single predefined integrator */
	{
		/* just one */
		fFields.Dimension(1);
		fFields = NULL;
		nControllerT* controller = fFEManager.nController(0);
		if (!controller) throw ExceptionT::kGeneralFail;
		
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
				field = new FieldT("displacement", NumSD(), *controller);

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
			{
				field = new FieldT("temperature", 1, *controller);				

				/* set labels */
				ArrayT<StringT> labels(1);
				labels[0] = "T";
				field->SetLabels(labels);
				break;
			}
			case GlobalT::kPML:
			{
				cout << "\n NodeManagerT::EchoFields: don't know how to configure\n"
				     <<   "     fields for PML analysis code: " << GlobalT::kPML << endl;
				throw ExceptionT::kGeneralFail;
			}			
			default:
				cout << "\nFEManagerT::SetController: unknown controller type\n" << endl;
				throw ExceptionT::kBadInputValue;
		}
		
		/* check */
		if (!field) throw ExceptionT::kGeneralFail;
		
		/* just one group */
		field->SetGroup(0);
		
		/* dimension */
		field->Dimension(NumNodes());
		
		/* clear all equation numbers */
		field->Equations() = FieldT::kInit;

		/* coordinate update field */
		if (field->Name() == "displacement") {
			if (fCoordUpdate) {
				cout << "\n NodeManagerT::EchoFields: \"displacement\" field already set" << endl;
				throw ExceptionT::kBadInputValue;
			}
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

void NodeManagerT::EchoExternalNodes(ostream& out)
{
	/* get external nodes */
	fFEManager.IncomingNodes(fExNodes);
	
	int num_ex = fExNodes.Length();
	if (num_ex > 0)
	{
		/* echo nodes */
		if (fFEManager.PrintInput())
		{
			out << " External nodes:\n";
			fExNodes++;
			out << fExNodes.wrap(8) << '\n';
			fExNodes--;
		}

		/* set up fields */
		for (int i = 0; i < fFields.Length(); i++)
			fFields[i]->InitExternalEquations(fExNodes);
	}
}

void NodeManagerT::EchoForceBC(FieldT& field, ifstreamT& in, ostream& out)
{
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
	{
		cout << "\n NodeManagerT::EchoForceBC: expecting numFBC_sets > 0: " 
		     << numFBCsets << endl;
		throw ExceptionT::kBadInputValue;
	}

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
		if (cards[i].Node() < 0 || cards[i].Node() >= NumNodes()) {
			cout << "\n NodeManagerT::EchoForceBC: node number is out of range: " 
			     << cards[i].Node() + 1 << endl;
			throw ExceptionT::kOutOfRange;
		}
}

void NodeManagerT::EchoHistoryNodes(ifstreamT& in, ostream &out)
{
	out << "\n N o d a l   H i s t o r i e s :\n\n";

	/* model database */
	ModelManagerT* model = fFEManager.ModelManager();

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
			if (!fExNodes.HasValue(node_set[j]))
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
		if (!model->RegisterElementGroup(ID, set, GeometryT::kPoint, true)) {
			cout << "\n NodeManagerT::EchoHistoryNodes: error initializing node set "
			     << node_set_ID[i] << " as model element ID " << ID << endl;
			throw ExceptionT::kBadInputValue;
		}
		fHistoryNodeSetIDs[i] = ID;
	}
}

#if 0
/* apply kinematic boundary conditions */
void NodeManagerT::ApplyKinematicBC(void)
{
	/* apply BC's */
	for (int i = 0; i < fKBC.Length(); i++)
		fnController->ConsistentKBC(fKBC[i]);
		
	/* KBC controllers */
	for (int j = 0; j < fKBCControllers.Length(); j++)
	{
		/* boundary condition cards generated by the controller */
		const ArrayT<KBC_CardT>& cards = fKBCControllers[j]->KBC_Cards();

		/* apply BC's */
		for (int i = 0; i < cards.Length(); i++)
			fnController->ConsistentKBC(cards[i]);
	}
}
#endif

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
		if (rowlabels->Length() != data.MajorDim()) throw ExceptionT::kSizeMismatch;
	
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
		case KBC_ControllerT::kSymmetricNodes:
		{
			SymmetricNodesT* kbc = new SymmetricNodesT(*this, field);
			return kbc;
		}
		case KBC_ControllerT::kPeriodicNodes:
		{
			PeriodicNodesT* kbc = new PeriodicNodesT(*this, field);
			return kbc;
		}
		default:
			cout << "\n NodeManagerT::NewKBC_Controller: KBC controller code "
			     << code <<   "     is not supported" << endl;
			throw ExceptionT::kBadInputValue;
	}
	return NULL;
}

FBC_ControllerT* NodeManagerT::NewFBC_Controller(FieldT& field, int code)
{
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
			fbc = new PenaltyWallT(fFEManager, field.Group(), eqnos, coords, velocity);
			break;

		case FBC_ControllerT::kAugLagWall:
			fbc = new AugLagWallT(fFEManager, this, field, coords);
			break;

		case FBC_ControllerT::kPenaltySphere:	
			fbc = new PenaltySphereT(fFEManager, field.Group(), eqnos, coords, velocity);
			break;

		case FBC_ControllerT::kAugLagSphere:	
			fbc = new AugLagSphereT(fFEManager, this, field, coords);
			break;

		case FBC_ControllerT::kMFPenaltySphere:	
			fbc = new MFPenaltySphereT(fFEManager, field.Group(), eqnos, coords, velocity);
			break;

		default:
			cout << "\n NodeManagerT::NewFBC_Controller: FBC controller code "
			     << code <<   "     is not supported" << endl;
			throw ExceptionT::kBadInputValue;
	}
	
	/* set time integrator */
	if (fbc) {
		const nControllerT& n_cont = field.nController();
		const eControllerT* e_cont = dynamic_cast<const eControllerT*>(&n_cont);
		fbc->SetController(e_cont);
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

	/* controller list */
	pArrayT<KBC_ControllerT*>& controllers = field.KBC_Controllers();

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
	controllers.Dimension(numKBC);
	for (int i = 0; i < numKBC; i++)
	{
		int num, type;
		in2 >> num >> type; num--;
		
		/* construct */
		controllers[num] = NewKBC_Controller(field, type);
		if (!controllers[num]) throw ExceptionT::kOutOfMemory;
		
		/* initialize */
		controllers[num]->Initialize(in2);
	}

	/* echo parameters */
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

	/* controller list */
	pArrayT<FBC_ControllerT*>& controllers = field.FBC_Controllers();

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

	controllers.Dimension(numFBCcont);
	for (int i = 0; i < numFBCcont; i++)
	{
		int num, type;
		in2 >> num >> type; num--;
		
		/* construct */
		controllers[num] = NewFBC_Controller(field, type);
		if (!controllers[num]) throw ExceptionT::kOutOfMemory;
		
		/* echo data */
		controllers[num]->EchoData(in2, out);

		/* initialize */
		controllers[num]->Initialize();
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
			dArray2DT n_values(node_set.Length(), 2*ndof);
			dArray2DT e_values;
		
			/* nodes in set */		
			for (int j = 0; j < node_set.Length(); j++) 
			{
				/* node map resolves processor-local node number to global
				 * node number. No map means serial execution */
				int node = node_set[j];
						
				/* compute reaction force */
				fFEManager.InternalForceOnNode(field, node, force);

				//TEMP - for now just write the displacements and the forces
				//       needs to be generalized for an arbitrary number of DOF
				//       and time derivatives.
				for (int k = 0; k < ndof; k++)
				{
					n_values(j, k) = field[0](node, k); /* "displacement" */
					n_values(j, k+ndof) = force[k];
				}
			}
			
			/* send for output */
			fFEManager.WriteOutput(ID, n_values, e_values);
		}
	}
}
