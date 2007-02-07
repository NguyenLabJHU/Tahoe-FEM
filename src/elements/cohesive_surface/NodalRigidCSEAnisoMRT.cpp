/* $Id: NodalRigidCSEAnisoMRT.cpp,v 1.3 2007-02-07 02:10:06 skyu Exp $ */
#include "NodalRigidCSEAnisoMRT.h"

#include "XDOF_ManagerT.h"
#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "SurfaceShapeT.h"
#include "InverseMapT.h"
#include "MR_NodalRP2DT.h"
#include "OutputSetT.h"

using namespace Tahoe;

//TEMP - debugging flag
#undef DEBUG

/* constructor */
NodalRigidCSEAnisoMRT::NodalRigidCSEAnisoMRT(const ElementSupportT& support):
	CSEAnisoT(support),
	fr(0),
	fCZRelation(NULL),
	fCurrPair(-1)
{
	SetName("nodal_rigid_anisotropic_CSE_MR");
}

/* destructor */
NodalRigidCSEAnisoMRT::~NodalRigidCSEAnisoMRT(void)
{

}

/* append element equations numbers to the list */
void NodalRigidCSEAnisoMRT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	const char caller[] = "NodalRigidCSEAnisoMRT::Equations";

	/* inherited */
	CSEAnisoT::Equations(eq_1, eq_2);
	
	/* collect the displacement equations */
	fCZNodePairDispEqnos.Dimension(fCZNodePairs.MajorDim(), 2*NumDOF());
	Field().SetLocalEqnos(fCZNodePairs, fCZNodePairDispEqnos);
	eq_1.AppendUnique(&fCZNodePairDispEqnos);

	/* source for XDOF information */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();

	/* resize the work space */
	int num_constraints = fConstraintXDOFTags.Length();
	fXDOFEqnos_man.SetMajorDimension(num_constraints, false);

	/* get the equations associate with the XDOF */
	const iArray2DT& xdof_eqnos = xdof_manager.XDOF_Eqnos(this, 0);
	if (xdof_eqnos.Length() != fConstraintXDOFTags.Length())
		ExceptionT::GeneralFail(caller, "expecting %d xdof equations not %d",
			fConstraintXDOFTags.Length(), xdof_eqnos.Length());
	
	/* displacement equations */
	const iArray2DT& disp_equations = Field().Equations();
	
	/* collect equations for active constraints */
	int constraint_dex = 0;
	int nprs = fConstraintStatus.MajorDim();
	int ndof = NumDOF();
	ArrayT<char> status;
	iArrayT constraint_eqnos;
	for (int i = 0; i < nprs; i++)
	{
		/* node pair */
		int n0 = fCZNodePairs(i,0);
		int n1 = fCZNodePairs(i,1);
	
		fConstraintStatus.RowAlias(i, status);
		for (int j = 0; j < ndof; j++)
		{
			/* collect equations */
			if (status[j] == kActive)
			{
				fXDOFEqnos.RowAlias(constraint_dex, constraint_eqnos);

				/* displacement equations */
				constraint_eqnos[0] = disp_equations(n0,j);
				constraint_eqnos[1] = disp_equations(n1,j);

				/* constraint equation */
				constraint_eqnos[2] = xdof_eqnos[constraint_dex];
				
				/* next */
				constraint_dex++;
			}
		}
	}

	/* add to list */
	eq_1.Append(&fXDOFEqnos);	
}

/* close current time increment */
void NodalRigidCSEAnisoMRT::CloseStep(void)
{
	/* inherited */
	CSEAnisoT::CloseStep();

	/* update constraint history */
	fConstraintStatus_n = fConstraintStatus;
	fStateVariables_n = fStateVariables;
	fConstraints_n = ElementSupport().XDOF_Manager().XDOF(this, 0);

	/* reset the 'last' state */
	fConstraintStatus_last = fConstraintStatus_n;
	fConstraints_last = fConstraints_n;
}

/* determine the number of constraints needed */
void NodalRigidCSEAnisoMRT::SetDOFTags(void)
{
	/* count number of active constraints */
	int count = 0;
	char* pa = fConstraintStatus.Pointer();
	int  len = fConstraintStatus.Length();
	for (int i = 0; i < len; i++)
		if (*pa++ == kActive)
			count++;

	/* resize tags array */
	fConstraintXDOFTags_man.SetLength(count, false);
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int NodalRigidCSEAnisoMRT::Reconfigure(void)
{
	/* store snapshot */
	fConstraintStatus_last = fConstraintStatus;
	fConstraints_last = ElementSupport().XDOF_Manager().XDOF(this, 0);

	/* run through constraints and check response */
	dArrayT state;
	int ndof = NumDOF();
	ArrayT<bool> is_rigid(ndof);
	dArrayT active_flags;
	dArrayT dummy;
	int num_changed = 0;
	for (int j = 0; j < fCZNodePairs.MajorDim(); j++)
	{
		/* query CZ response */
		fStateVariables.RowAlias(j, state);		
		fCZRelation->RigidQ(state, is_rigid);

		/* get the active flags */
		fCZRelation->GetActiveFlags(state, active_flags);
			
		/* check if relations have reached failure */
		SurfacePotentialT::StatusT status = fCZRelation->Status(dummy, state);

		/* loop over dof */
		for (int i = 0; i < ndof; i++)
		{
			/* check to see if either of the nodes is unprescribed */
			bool no_free_disp = false;
			if (fCZNodePairDispEqnos(j,i) < 1 || fCZNodePairDispEqnos(j,ndof+i) < 1)
				no_free_disp = true;

			char& curr_state = fConstraintStatus(j,i);
			if (status == SurfacePotentialT::Failed && curr_state != kFailed)
			{
				curr_state = kFailed;
				num_changed++;
				active_flags[i] = 1.0; /* evolution active */
			}
			else if (no_free_disp) /* neither node has a free displacement */
			{
				/* remove the constraint */
				if (curr_state == kActive)
				{
					curr_state = kUnneeded;
					num_changed++;
					active_flags[i] = 0.0; /* don't solve for traction */
				}
			}
			else if (status == SurfacePotentialT::Precritical) /* not enough stress to initiate */
			{
				/* activate constraint */
				if (curr_state == kFree) {
					curr_state = kActive;
					num_changed++;
					active_flags[i] = 0.0; /* evolution not active */
				}
			}
			else if (curr_state == kActive && is_rigid[i] == false)
			{
				/* release constraint */
				curr_state = kFree;
				num_changed++;
				active_flags[i] = 1.0; /* evolution active */
			}
			else if (curr_state == kFree && is_rigid[i] == true)
			{
				/* enable constraint */
				curr_state = kActive;
				num_changed++;
				active_flags[i] = 0.0; /* evolution not active */
			}
		}
	}

	/* signal changes */
	if (num_changed > 0)
		return 1;
	else
		return 0;
}

/* restore any state data to the previous converged state. */
void NodalRigidCSEAnisoMRT::ResetState(void)
{
	/* inherited */
	CSEAnisoT::ResetStep();

	/* reset the status */
	fStateVariables = fStateVariables_n;

	/* reset the 'last' state */
	fConstraintStatus_last = fConstraintStatus_n;
	fConstraints_last = fConstraints_n;
}

/* collecting element connectivities for the field */
void NodalRigidCSEAnisoMRT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1, 
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	CSEAnisoT::ConnectsU(connects_1, connects_2);

	/* add connectivities with constraints */
	connects_1.AppendUnique(&fXDOFConnectivities);
}

iArrayT& NodalRigidCSEAnisoMRT::DOFTags(int tag_set)
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("NodalRigidCSEAnisoMRT::DOFTags", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif
	return fConstraintXDOFTags;
}

/* return the contact elements */
const iArray2DT& NodalRigidCSEAnisoMRT::DOFConnects(int tag_set) const
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("NodalRigidCSEAnisoMRT::DOFConnects", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif
	return fXDOFConnectivities;
}

/* restore the DOF values to the last converged solution */
void NodalRigidCSEAnisoMRT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("NodalRigidCSEAnisoMRT::ResetDOF", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif

	dArrayT state, traction;
	int ndof = NumDOF();
	int pair = 0;
	int dof = 0;
	int curr_dof_count = 0;
	int last_dof_count = 0;
	int num_status = fConstraintStatus.Length();
	for (int i = 0; i < num_status; i++)
	{
		char curr_status = fConstraintStatus[i];
		char last_status = fConstraintStatus_last[i];
	
		if (curr_status == kActive) {
			if (last_status == kActive)
				DOF[curr_dof_count] = fConstraints_last[last_dof_count];
			else /* initialize constraint from traction */
			{
				/* traction vector */
				fStateVariables.RowAlias(pair, state);
				fCZRelation->GetTraction(state, traction);
			
				/* force */
				DOF[curr_dof_count] = traction[dof]*fCZNodeAreas[pair];
			}
			curr_dof_count++;
		}
		
		if (last_status == kActive)
			last_dof_count++;
			
		if (++dof == ndof) {
			pair++;
			dof = 0;
		}
	}
}

/* generate nodal connectivities */
void NodalRigidCSEAnisoMRT::GenerateElementData(void)
{
	/* resize the work space */
	int num_constraints = fConstraintXDOFTags.Length();
	fXDOFConnectivities_man.SetMajorDimension(num_constraints, false);

	/* collect nodes and DOF tags */
	int num_element_constraints = fConstraintStatus.MinorDim();
	int tag = 0;
	int num_pairs = fCZNodePairs.MajorDim();
	iArrayT row, pair;
	for (int i = 0; i < num_pairs; i++)
	{
		fCZNodePairs.RowAlias(i, pair);
		char* constraint_status = fConstraintStatus(i);
		for (int j = 0; j < num_element_constraints; j++)
			if (*constraint_status++ == kActive)
			{
				/* collect connectivities */
				fXDOFConnectivities.RowAlias(tag, row);
				row[0] = pair[0];
				row[1] = pair[1];
				row[2] = fConstraintXDOFTags[tag];
		
				/* next */
				tag++;
			}
	}

	/* echo connectivities */
	if (ElementSupport().PrintInput()) {
		ofstreamT& out = ElementSupport().Output();
		out << "\n NodalRigidCSEAnisoMRT::GenerateElementData: constraint connectivities\n" << '\n';
		fXDOFConnectivities++;
		fXDOFConnectivities.WriteNumbered(out);
		fXDOFConnectivities--;
	}
}

void NodalRigidCSEAnisoMRT::RegisterOutput(void)
{
	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	
	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

	/* set output specifier */
	OutputSetT output_set(GeometryT::kPoint, fCZNodePairPoints, n_labels);

	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}
	
void NodalRigidCSEAnisoMRT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	/* dimensions */
	int  nsd = NumSD();
	int ndof = NumDOF();
	int  nen = NumElementNodes();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* allocate results space */
	e_values.Dimension(0,0);
	
	/* all values per node */
	RaggedArray2DT<double> node_all;
	node_all.Configure(n_codes);
	dArrayT coords, disp, jump, T, matdat;	
	node_all.RowAlias(0, coords);
	node_all.RowAlias(1, disp);
	node_all.RowAlias(2, jump);
	node_all.RowAlias(3, T);
	node_all.RowAlias(4, matdat);

	/* reference coordinates */
	const dArray2DT& ref_coords = ElementSupport().InitialCoordinates();

	/* displacement field */
	const FieldT& field = Field();
	const dArray2DT& disp_field = field[0];

	/* constraint forces */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();
	const dArray2DT& constraints = xdof_manager.XDOF(this, 0);

	/* for assembly */
	iArrayT node(1);
	dArray2DT nodal_values(1, node_all.Length(), node_all.Pointer());

	/* loop over pairs */
	int constraint_dex = 0;
	dArrayT d0, d1, state, diff_d0d1(NumDOF()), sigma;
	ArrayT<char> constraint_state;
	for (int i = 0; i < fCZNodePairs.MajorDim(); i++)
	{
		/* pair nodes */
		int n0 = fCZNodePairs(i,0);
		int n1 = fCZNodePairs(i,1);
		
		/* coordinate transformation */
		// double Q = fCZDirection[i];
		double j1 = fCurrShapes->Jacobian(Q);

	
		/* coordinates - same for both nodes in pair */
		if (n_codes[NodalCoord])
			ref_coords.RowAlias(n0, coords);
	
		/* displacements and jump */
		disp_field.RowAlias(n0, d0);
		disp_field.RowAlias(n1, d1);
		diff_d0d1.DiffOf(d1,d0);
		// diff_d0d1[0] *= Q;
		// diff_d0d1[1] *= Q;
		Q.MultTx(diff_d0d1, diff_d0d1);

		if (n_codes[NodalDispJump])
			jump = diff_d0d1;

		/* state variables */
		fStateVariables.RowAlias(i,state);
		fConstraintStatus.RowAlias(i, constraint_state);
		
		/* traction */
		if (n_codes[NodalTraction]) {
		
			/* traction from cohesive relations */
			T = fCZRelation->Traction(diff_d0d1, state, sigma, false);
		
			/* overwrite with active constraints */
			for (int j = 0; j < ndof; j++)
				if (constraint_state[j] == kActive) /* rigid -> traction from constraints */
				{
					T[j] = constraints[constraint_dex]/fCZNodeAreas[i];
					constraint_dex++;
				}
		}

		/* material output data */
		if (n_codes[MaterialData])
			fCZRelation->ComputeOutput(diff_d0d1, state, matdat);
		
		/* first node */
		if (n_codes[NodalDisp]) 
			disp = d0;
		node[0] = n0;
		ElementSupport().AssembleAverage(node, nodal_values);
	
		/* second node */
		if (n_codes[NodalDisp]) 
			disp = d1;
		node[0] = n1;
		ElementSupport().AssembleAverage(node, nodal_values);
	}

	/* get ordered, nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
}

/* write restart data to the output stream */
void NodalRigidCSEAnisoMRT::WriteRestart(ostream& out) const
{
	/* inherited */
	CSEAnisoT::WriteRestart(out);

	/* last constraint status - fConstraintStatus_n is written/read instead
	 * of fConstraintStatus because upon restart fConstraintStatus associated
	 * with the reference configuration will be compared against fConstraintStatus_last
	 * to see whether the equation system needs to be reset, and since WriteRestart
	 * is called after CloseStep, fConstraintStatus_n = fConstraintStatus for
	 * the current state. */
	for (int i = 0; i < fConstraintStatus_n.Length(); i++)
		out << fConstraintStatus_n[i];
	out << '\n';

	/* previous value of the constraints */
	dArrayT tmp;
	tmp.Alias(fConstraints_n);
	out << tmp.Length() << '\n' << tmp.wrap_tight(5) << '\n';
	
	/* state variables */
	out << fStateVariables.wrap_tight(5) << '\n';
}

/* read restart data to the output stream */
void NodalRigidCSEAnisoMRT::ReadRestart(istream& in)
{
	/* inherited */
	CSEAnisoT::ReadRestart(in);

	/* constraint status */
	for (int i = 0; i < fConstraintStatus_n.Length(); i++)
		in >> fConstraintStatus_n[i];

	/* previous value of the constraints */
	int len;
	in >> len;
	fConstraints_n.Dimension(len);
	dArrayT tmp;
	tmp.Alias(fConstraints_n);
	in >> tmp;
	
	/* reset 'last' state */
	fConstraints_last = fConstraints_n;
	fConstraintStatus_last = fConstraintStatus_n;

	/* state variables */
	in >> fStateVariables;
	fStateVariables_n = fStateVariables;
}

/* describe the parameters needed by the interface */
void NodalRigidCSEAnisoMRT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	CSEAnisoT::DefineParameters(list);

	/* regularization */
	ParameterT regularization(fr, "regularization");
	regularization.SetDefault(fr);
	regularization.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(regularization);
}

/* accept parameter list */
void NodalRigidCSEAnisoMRT::TakeParameterList(const ParameterListT& list)
{
	/* regularization */
	fr = list.GetParameter("regularization");

	/* inherited */
	CSEAnisoT::TakeParameterList(list);

	const char caller[] = "NodalRigidCSEAnisoMRT::Initialize";

	/* inherited */
	fElementCards.Current(0); /* reset element counter even though it's not used */

	/* local node numbers on each facet */
	const iArray2DT& nodes_on_facets = fShapes->NodesOnFacets();
	int nfn = nodes_on_facets.MinorDim();
	iArrayT face1, face2;
	nodes_on_facets.RowAlias(0, face1);
	nodes_on_facets.RowAlias(1, face2);

	/* count cohesive zone nodes */
	ArrayT<char> cz_node(ElementSupport().NumNodes());
	cz_node = 0;
	int count = 0;
	int nel = NumElements();
	for (int i = 0; i < nel; i++) 
	{
		const ElementCardT& card = ElementCard(i);
		const iArrayT& nodes = card.NodesU();
		for (int j = 0; j < nfn; j++) 
		{
			int nd1 = nodes[face1[j]];
			int nd2 = nodes[face2[j]];
			if (nd1 != nd2 && cz_node[nd1] == 0) /* first time and no self-self constraints */
			{
				cz_node[nd1] = 1;
				count++;
			}
		}
	}
	fCZNodePairs.Dimension(count,2);
	fCZNodePairPoints.Alias(fCZNodePairs.Length(), 1, fCZNodePairs.Pointer());

	/* collect pairs */
	cz_node = 0;
	count = 0;
	for (int i = 0; i < nel; i++) 
	{
		const ElementCardT& card = ElementCard(i);
		const iArrayT& nodes = card.NodesU();
		for (int j = 0; j < nfn; j++) 
		{
			int nd_1 = nodes[face1[j]];
			int nd_2 = nodes[face2[j]];
			if (nd_1 != nd_2 && cz_node[nd_1] == 0) /* first time */
			{
				cz_node[nd_1] = 1;
				fCZNodePairs(count, 0) = nd_1;
				fCZNodePairs(count, 1) = nd_2;
				count++;
			}
		}
	}

	/* construct an inverse map to local of first node of each pair */
	iArrayT first_node(fCZNodePairs.MajorDim());
	fCZNodePairs.ColumnCopy(0, first_node);
	InverseMapT global_to_pair_map;
	global_to_pair_map.SetOutOfRange(InverseMapT::MinusOne);
	global_to_pair_map.SetMap(first_node);

	/* compute the triburary area of each node */
	fCZNodeAreas.Dimension(fCZNodePairs.MajorDim());
	fCZNodeAreas = 0.0;
	fCZDirection.Dimension(fCZNodePairs.MajorDim());
	fCZDirection = 0.0;

	/* allocate work space */
	// dMatrixT Q(NumSD());
	Q.Dimension(NumSD());;
	int nee = NumElementNodes()*NumDOF();
	dQ.Dimension(NumSD());

        for (int k = 0; k < NumSD(); k++)
		dQ[k].Dimension(NumSD(), nee);

	dArrayT area(NumElementNodes()), Na(NumElementNodes());
	for (int i = 0; i < nel; i++) 
	{
		/* element information */
		const ElementCardT& card = ElementCard(i);
		const iArrayT& nodes_X = card.NodesX();
	
		/* get ref geometry (1st facet only) */
		fNodes1.Collect(face1, nodes_X);
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords);

		/* integrate */
		area = 0.0;
		fShapes->TopIP();
		while (fShapes->NextIP())
		{  
			/* integration weights */
			double w = fShapes->IPWeight();	
			// double j0 = fShapes->Jacobian(Q);
			double j0 = fShapes->Jacobian();
			double j1 = fCurrShapes->Jacobian(Q, dQ);

			//TEMP - require normal to be in the x2 direction:
			//if (fabs(fabs(Q(1,1)) - 1.0) > kSmall)
			//	ExceptionT::GeneralFail(caller, "surface normals must be in the x2 direction");

			/* get shape functions */
			fShapes->Shapes(Na);
			
			/* integrate */
			area.AddScaled(w*j0, Na);
		}
		
		/* accumulate nodal values */
		for (int j = 0; j < area.Length(); j++)
		{
			int pair = global_to_pair_map.Map(nodes_X[j]);
			if (pair != -1) {
			
				/* nodal area */
				fCZNodeAreas[pair] += area[j];
				
				/* direction from the last ip */
			//	if (fabs(fCZDirection[pair]) < 0.5)
			//		fCZDirection[pair] = Q(1,1);
			}
		}
	}

	/* flags array */
	fConstraintStatus.Dimension(fCZNodePairs.MajorDim(), NumDOF());
	fConstraintStatus = kActive;
	fConstraintStatus_n = fConstraintStatus;
	fConstraintStatus_last = fConstraintStatus;

	/* equations in every constrained pair */
	int neq = 2 + 1; /* 2 nodes and the constraint */

	/* dynamic work space managers */
	fConstraintXDOFTags_man.SetWard(0, fConstraintXDOFTags);
	fXDOFConnectivities_man.SetWard(0, fXDOFConnectivities, neq);
	fXDOFEqnos_man.SetWard(0, fXDOFEqnos, neq);

	/* register with node manager - sets initial fContactDOFtags */
	iArrayT xdof_tags(1);
	xdof_tags = 1;
	ElementSupport().XDOF_Manager().XDOF_Register(this, xdof_tags);

	/* allocate state variables */
	if (fSurfPots.Length() != 1) ExceptionT::BadInputValue(caller, "only 1 potential supported: %d", fSurfPots.Length());
	SurfacePotentialT* surf_pot = fSurfPots[0];
	fCZRelation = dynamic_cast<MR_NodalRP2DT*>(surf_pot);
	if (!fCZRelation) ExceptionT::BadInputValue(caller, "cohesive relation must be MR_NodalRP2DT");
	int num_state = fCZRelation->NumStateVariables(); 
	fStateVariables.Dimension(fCZNodePairs.MajorDim(), num_state);

	/* initialize state variable space */
	dArrayT state;
	dArrayT active_flags;
	for (int i = 0; i < fCZNodePairs.MajorDim(); i++)
	{
		fStateVariables.RowAlias(i, state);
		fCZRelation->InitStateVariables(state);
		
		/* mark all evolution equations as inactive since all the constraints are active */
		fCZRelation->GetActiveFlags(state, active_flags);
		active_flags = 0.0;
	}

	/* set history */
	fStateVariables_n = fStateVariables;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* force vector */
void NodalRigidCSEAnisoMRT::RHSDriver(void)
{
	const char caller[] = "NodalRigidCSEAnisoMRT::RHSDriver";

	/* time-integration parameters */
	double constKd = 1.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* source for XDOF information */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();

	/* current values of the constraints */
	const dArray2DT& constraints = xdof_manager.XDOF(this, 0);
	if (constraints.Length() != fConstraintXDOFTags.Length())
		ExceptionT::GeneralFail(caller, "expecting %d constraints not %d",
			fConstraintXDOFTags.Length(), constraints.Length());

	/* displacements */
	const FieldT& field = Field();
	const dArray2DT& disp = field(0,0);
	const dArray2DT& disp_last = field(-1,0);

	/* set state to start of current step */
	fStateVariables = fStateVariables_n;

	/* compute nodal forces */
	int ndof = NumDOF();
	dArrayT jump(ndof), jump_last(ndof);
	dArrayT rhs(fXDOFEqnos.MinorDim()), rhs_disp(2*ndof);
	int constraint_dex = 0;
	dArrayT state, traction, ftraction;
	dArrayT dummy;
	ArrayT<char> constraint_status;
	iArrayT eqnos;
	dArrayT h(ndof), f(ndof), l(ndof);
	dMatrixT QT(NumSD());
	for (int j = 0; j < fCZNodePairs.MajorDim(); j++)
	{
		/* node pair */
		fCurrPair = j;
		int n0 = fCZNodePairs(j,0);
		int n1 = fCZNodePairs(j,1);
	
		/* element state */
		fConstraintStatus.RowAlias(j, constraint_status);
		fStateVariables.RowAlias(j, state);
		fCZRelation->GetTraction(state, traction);

		/* compute displacement jumps */
		jump.DiffOf(disp(n1), disp(n0));
		jump_last.DiffOf(disp_last(n1), disp_last(n0));

		/* coordinate transformation */
		// double Q = fCZDirection[j];
		double j1 = fCurrShapes->Jacobian(Q);

		// jump[0] *= Q;
		// jump[1] *= Q;
		// jump_last[0] *= Q;
		// jump_last[1] *= Q;
		Q.MultTx(jump, jump);
		Q.MultTx(jump_last, jump_last);

		/* loop over dof */
		int num_constrained = 0;
		int num_failed = 0;
		int max_num_constrained = ndof;

		for (int i = 0; i < ndof; i++)
			if (constraint_status[i] == kActive) /* enforce constraint */
			{
				/* constraint and multiplier */
				// double h = jump[i] - jump_last[i];
				// double f = constraints[constraint_dex];
				// double l = f + fr*h;
				h[i] = jump[i] - jump_last[i];
				f[i] = constraints[constraint_dex];
				l[i] = f[i] + fr*h[i];

				/* update traction in the state variable array */
				// traction[i] = f/fCZNodeAreas[j];
				traction[i] = f[i]/fCZNodeAreas[j];

				/* residual of the displacement equations */
				// rhs[0] = l*Q; /* -1*-l */
				// rhs[1] =-l*Q; /* +1*-l */
				rhs[0] = l[i]*Q(i,i);
				rhs[1] =-l[i]*Q(i,i);

				/* residual of the constraint equation */
				// rhs[2] =-h;
				rhs[2] =-h[i]*Q(i,i);

				if (i == 0 && constraint_status[1] == kActive)
				{
					h[1] = jump[1] - jump_last[1];
					f[1] = constraints[constraint_dex+1];
					l[1] = f[1] + fr*h[1];

					rhs[0] += l[1]*Q(1,0);
					rhs[1] +=-l[1]*Q(1,0);
					rhs[2] +=-h[1]*Q(1,0);
				}

				if (i == 1 && constraint_status[0] == kActive)
				{
					h[0] = jump[0] - jump_last[0];
					f[0] = constraints[constraint_dex-1];
					l[0] = f[0] + fr*h[0];

					rhs[0] += l[0]*Q(0,1);
					rhs[1] +=-l[0]*Q(0,1);
					rhs[2] +=-h[0]*Q(0,1);
				}

				/* get equation numbers */
				fXDOFEqnos.RowAlias(constraint_dex, eqnos);
	
				/* assemble */
				ElementSupport().AssembleRHS(Group(), rhs, eqnos);

				/* next constraint */
				constraint_dex++;
				num_constrained++;
			}
			else if (constraint_status[i] == kFailed)
				num_failed++;
			else if (constraint_status[i] == kUnneeded)
				max_num_constrained--;

		/* at least one active cz relation */
		if (num_failed == 0 && num_constrained < max_num_constrained)
		{
			/* call traction */
			const dArrayT& traction = fCZRelation->Traction(jump, state, dummy, true);

			Q.Multx(traction, ftraction);

			/* expand */
			fShapes->Grad_d().MultTx(ftraction, ftraction); // need to be checked

			/* loop over directions */
			for (int i = 0; i < ndof; i++)
				if (constraint_status[i] == kActive)
				{
					/* no cz contribution */
					rhs_disp[i] = 0.0;
					rhs_disp[i+ndof] = 0.0;
				}
				else /* contributon from cz tractions */
				{
					/* equivalent nodal force */
					// double f =-Q*traction[i]*fCZNodeAreas[j];
					double f1 =-ftraction[i]*fCZNodeAreas[j];

					/* equal and opposite */
					rhs_disp[i] = -f1;
					rhs_disp[i+ndof] = f1;
				}
#ifdef DEBUG
ElementSupport().Output() << "jmp = " << jump.no_wrap() << '\n';
ElementSupport().Output() << "trc = " << traction.no_wrap() << '\n';
ElementSupport().Output() << "rhs = " << rhs_disp.no_wrap() << endl;
#endif	
			/* get equations numbers */
			fCZNodePairDispEqnos.RowAlias(j, eqnos);
			
			/* assemble */
			ElementSupport().AssembleRHS(Group(), rhs_disp, eqnos);
		}
	}
}

/* tangent matrix */
void NodalRigidCSEAnisoMRT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	const char caller[] = "NodalRigidCSEAnisoMRT::LHSDriver";

	/* matrix format */
	dMatrixT::SymmetryFlagT format = (fRotate) ?
		dMatrixT::kWhole : dMatrixT::kUpperOnly;


	/* time-integration parameters */
	double constK = 1.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* jumps */
	int ndof = NumDOF();

	/* workspace */
	dArrayT vec_nee(fXDOFEqnos.MinorDim());
	vec_nee = 0.0;
	// ElementMatrixT lhs(fXDOFEqnos.MinorDim(), ElementMatrixT::kSymmetric);
	// ElementMatrixT lhs_disp(2*ndof, ElementMatrixT::kSymmetric);
	ElementMatrixT lhs(fXDOFEqnos.MinorDim(), ElementMatrixT::kNonSymmetric);
	ElementMatrixT lhs_disp(2*ndof, ElementMatrixT::kNonSymmetric);
	dMatrixT stiffness(ndof), djump_du(ndof, 2*ndof);
	djump_du = 0.0;

	/* displacements */
	const FieldT& field = Field();
	const dArray2DT& disp = field(0,0);

	/* constraint stiffness matrix is constant */
	// lhs(0,0) = fr;
	// lhs(1,0) =-fr;
	// lhs(0,1) =-fr;
	// lhs(1,1) = fr;
	// lhs(2,2) = 0;

	int constraint_dex = 0;
	ArrayT<char> constraint_status;
	iArrayT eqnos;
	dArrayT jump(ndof), state, sigma, ftraction;
	for (int j = 0; j < fCZNodePairs.MajorDim(); j++)
	{
		/* node pair */
		fCurrPair = j;
		int n0 = fCZNodePairs(j,0);
		int n1 = fCZNodePairs(j,1);

		/* pair state */
		fConstraintStatus.RowAlias(j, constraint_status);
		fStateVariables.RowAlias(j, state);

		/* coordinate transformation */
		// double Q = fCZDirection[j];
		double j1 = fCurrShapes->Jacobian(Q, dQ);

		/* loop over dof */
		int num_constrained = 0;
		int num_failed = 0;
		for (int i = 0; i < ndof; i++)
			if (constraint_status[i] == kActive)
			{
				/* get equation numbers */
				fXDOFEqnos.RowAlias(constraint_dex, eqnos);

				/* constraint gradient */
				// lhs(2,0) =-Q;
				// lhs(2,1) = Q;
				// lhs(0,2) =-Q;
				// lhs(1,2) = Q;
				lhs(0,0) = fr*Q(0,i);
				lhs(1,0) =-fr*Q(0,i);
				lhs(2,0) =-Q(0,i);
				lhs(0,1) = fr*Q(1,i);
				lhs(1,1) =-fr*Q(1,i);
				lhs(2,1) =-Q(1,i);
				lhs(0,2) = 1;
				lhs(1,2) =-1;
				lhs(2,2) = 0;

				/* assemble */
				ElementSupport().AssembleLHS(Group(), lhs, eqnos);
			
				/* next constraint */
				constraint_dex++;
				num_constrained++;
			}
			else if (constraint_status[i] == kFailed)
				num_failed++;

		/* at least one active cz relation */
		if (num_failed == 0 && num_constrained < ndof)
		{		
			/* loop over directions */
			djump_du = fShapes->Grad_d();   //need to be checked

			for (int i = 0; i < ndof; i++)
			{
				if (constraint_status[i] == kActive)
				{
					djump_du(i,i) = 0.0;
					djump_du(i,i+ndof) = 0.0;
				}
				else /* contributon from cz tractions */
				{
					// djump_du(i,i) = -Q;
					// djump_du(i,i+ndof) = Q;
				}
			}

			/* call stiffness */
			jump.DiffOf(disp(n1), disp(n0));
			// jump[0] *= Q;
			// jump[1] *= Q;
			Q.MultTx(jump, jump);

			/* stiffness in local frame */
			const dMatrixT& K = fCZRelation->Stiffness(jump, state, sigma);

			/* rotation */
			/* traction in local frame */
			const dArrayT& traction = fCZRelation->Traction(jump, state, sigma, false);

			/* 1st term */
			ftraction.SetToScaled(fCZNodeAreas[j], traction);
			CSEAnisoT::Q_ijk__u_j(dQ, ftraction, fnsd_nee_1);
			fNEEmat.MultATB(djump_du, fnsd_nee_1);
			lhs_disp += fNEEmat;

			/* 2nd term */
			stiffness.SetToScaled(fCZNodeAreas[j], K);
			fnsd_nee_1.MultATB(Q, djump_du);
			fnsd_nee_2.MultATB(stiffness, fnsd_nee_1);
			u_i__Q_ijk(jump, dQ, fnsd_nee_1);
			fNEEmat.MultATB(fnsd_nee_2, fnsd_nee_1);
			lhs_disp += fNEEmat;

			/* 3rd term */
			stiffness.MultQBQT(Q, K);
			// stiffness.SetToScaled(fCZNodeAreas[j], fCZRelation->Stiffness(jump, state, sigma));
			stiffness.SetToScaled(fCZNodeAreas[j], stiffness);

			/* compute element stiffness */
			// lhs_disp.MultQTBQ(djump_du, stiffness);
			lhs_disp.MultQTBQ(djump_du, stiffness, format, dMatrixT::kAccumulate);


//TEMP
#ifdef DEBUG
ofstreamT& out = ElementSupport().Output();
out << "K =\n" << lhs_disp  << endl;
#endif

			/* get equation numbers */
			fCZNodePairDispEqnos.RowAlias(j, eqnos);
	
			/* assemble */
			ElementSupport().AssembleLHS(Group(), lhs_disp, eqnos);
		}		
	}
}

/* write all current element information to the stream */
void NodalRigidCSEAnisoMRT::CurrElementInfo(ostream& out) const
{
	/* node pair */
	int n0 = fCZNodePairs(fCurrPair, 0);
	int n1 = fCZNodePairs(fCurrPair, 1);

	/* reference coordinates */
	const dArray2DT& ref_coords = ElementSupport().InitialCoordinates();
	
	/* pair state */
	ArrayT<char> status;
	fConstraintStatus.RowAlias(fCurrPair, status);
	dArrayT state;
	fStateVariables.RowAlias(fCurrPair, state);
	
	/* displacement field */
	const FieldT& field = Field();
	const dArray2DT& disp = field(0,0);
	dArrayT jump(NumDOF());
	jump.DiffOf(disp(n1), disp(n0));

	/* write info */
	out << "\nNodalRigidCSEAnisoMRT::CurrElementInfo\n";
	out << "pair number = " << fCurrPair+1 << '\n';
	out << "n_1 = " << n0+1 << " : ";
	ref_coords.PrintRow(n0, out);
	out << "n_2 = " << n1+1 << " : ";
	ref_coords.PrintRow(n1, out);
	for (int i = 0; i < status.Length(); i++)
		out << i+1 << " : " << ((status[i] == kActive) ? "active" : "free") << '\n';
	out << "state =\n" << state.wrap(5) << '\n';
	out << " jump = " << jump.no_wrap() << endl; 
}
