/* $Id: FEManagerT_bridging.cpp,v 1.3.2.15 2003-05-14 19:49:54 hspark Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "KBC_PrescribedT.h"
#include "KBC_CardT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "NLSolver.h"

#include "BridgingScaleT.h"
#include "ParticleT.h"
#include "dSPMatrixT.h"

using namespace Tahoe;

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT(input, output, comm),
	fBridgingIn(bridging_input),
	fParticle(NULL),
	fBridgingScale(NULL),
	fSolutionDriver(NULL)
{

}

/* send update of the solution to the NodeManagerT */
void FEManagerT_bridging::Update(int group, const dArrayT& update)
{
	/* accumulative */
	dArrayT& cumulative_update = fCumulativeUpdate[group];
	if (cumulative_update.Length() == update.Length())
		cumulative_update += update;
	
	/* inherited */
	FEManagerT::Update(group, update);
}

/* compute RHS-side, residual force vector and assemble to solver */
void FEManagerT_bridging::FormRHS(int group) const
{
	/* inherited */
	FEManagerT::FormRHS(group);

	/* assemble external contribution */
	const dArrayT* external_force = fExternalForce[group];
	if (external_force != NULL) {
		fSolvers[group]->UnlockRHS();
		fSolvers[group]->AssembleRHS(*external_force);
		fSolvers[group]->LockRHS();
	}
	
	/* assemble external contribution */
	const dArray2DT* external_force_2D = fExternalForce2D[group];
	if (external_force_2D != NULL) {
		fSolvers[group]->UnlockRHS();
		fSolvers[group]->AssembleRHS(*external_force_2D, fExternalForce2DEquations[group]);
		fSolvers[group]->LockRHS();		
	}
}

/* reset the cumulative update vector */
void FEManagerT_bridging::ResetCumulativeUpdate(int group)
{
	fCumulativeUpdate[group].Dimension(fNodeManager->NumEquations(group));
	fCumulativeUpdate[group] = 0.0;
}

/* (re-)set the equation number for the given group */
void FEManagerT_bridging::SetEquationSystem(int group)
{
	/* inherited */
	FEManagerT::SetEquationSystem(group);

	//NOTE: this is going to break if the equation numbers has changed since the force was set
	if (fExternalForce2D[group])
		ExceptionT::GeneralFail("FEManagerT_bridging::SetEquationSystem",
			"group %d has external force so equations cannot be reset", group+1);
}

/* set pointer to an external force vector */
void FEManagerT_bridging::SetExternalForce(const StringT& field, const dArray2DT& external_force, const iArrayT& activefenodes)
{
	const char caller[] = "FEManagerT_bridging::SetExternalForce";

	/* check */
	if (activefenodes.Length() != external_force.MajorDim()) 
		ExceptionT::SizeMismatch(caller);

	/* get the field */
	const FieldT* thefield = fNodeManager->Field(field);
	if (!thefield) ExceptionT::GeneralFail(caller);

	/* store pointers */
	int group = thefield->Group();
	fExternalForce2D[group] = &external_force;
	fExternalForce2DNodes[group] = &activefenodes;
	
	/* collect equation numbers */
	iArray2DT& eqnos = fExternalForce2DEquations[group];
	eqnos.Dimension(fNodeManager->NumNodes(), thefield->NumDOF());
	thefield->SetLocalEqnos(activefenodes, eqnos);	// crashes here if not nodes and atoms everywhere
}

/* initialize the ghost node information */
void FEManagerT_bridging::InitGhostNodes(void)
{
	const char caller[] = "FEManagerT_bridging::InitGhostNodes";

	/* collect ghost nodes */
	if (fBridgingIn.is_open()) {
		ArrayT<StringT> id_list;
		fModelManager->NodeSetList(fBridgingIn, id_list);
		fModelManager->ManyNodeSets(id_list, fGhostNodes);
	}

	/* assume atomistic field is "displacement" */
	StringT field = "displacement";
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* create controller to driver solution */
	if (!fSolutionDriver) {
	
		/* construct new contoller */
		fSolutionDriver = new KBC_PrescribedT(*fNodeManager);
	
		/* add to field */
		the_field->AddKBCController(fSolutionDriver);
	}

	/* generate KBC cards - all degrees of freedom */
	int ndof = the_field->NumDOF();
	ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
	KBC_cards.Dimension(fGhostNodes.Length()*ndof);
	int dex = 0;
	for (int j = 0; j < ndof; j++)
		for (int i = 0; i < fGhostNodes.Length(); i++)
			KBC_cards[dex++].SetValues(fGhostNodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	
	/* reset neighbor lists */
	ParticleT& particle = Particle();
	particle.SetSkipParticles(fGhostNodes);
	particle.SetConfiguration();
	
	/* reset the group equations numbers */
	SetEquationSystem(the_field->Group());

	/* echo ghost nodes */
	if (fPrintInput) {
		fGhostNodes++;
		fMainOut << "\n Ghost nodes:\n";
		fMainOut << fGhostNodes.wrap(5) << '\n';
		fGhostNodes--;
	}

	/* mark nodes as ghost */
	fNonGhostNodes.Dimension(fModelManager->NumNodes() - fGhostNodes.Length());
	iArrayT is_ghost(fModelManager->NumNodes());
	is_ghost = 0;
	for (int i = 0; i < fGhostNodes.Length(); i++)
		is_ghost[fGhostNodes[i]] = 1;

	/* check for uniqueness */
	int ng = is_ghost.Count(1);
	if (ng != fGhostNodes.Length())
		ExceptionT::GeneralFail(caller, "list of ghost nodes contains %d duplicates",
			fGhostNodes.Length() - ng);

	/* collect non-ghost nodes */
	dex = 0;
	for (int i = 0; i < is_ghost.Length(); i++)
		if (is_ghost[i] == 0)
			fNonGhostNodes[dex++] = i;
}

/* compute the ghost-nonghost part of the stiffness matrix */
void FEManagerT_bridging::Form_G_NG_Stiffness(const StringT& field, dSPMatrixT& K_G_NG)
{
	const char caller[] = "FEManagerT_bridging::Form_G_NG_Stiffness";

	/* get the field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* redimension if needed */
	int row_eq = the_field->NumEquations();
	int col_eq = fGhostNodes.Length()*the_field->NumDOF();
	if (K_G_NG.Rows() != row_eq || K_G_NG.Cols() != col_eq) {
		K_G_NG.Dimension(row_eq, col_eq, 0);
	}

	/* dimension pseudo equations array and map */
	if (fGhostNodesEquations.MajorDim() != fGhostNodes.Length() ||
	    fGhostNodesEquations.MinorDim() != the_field->NumDOF()) {
		fGhostNodesEquations.Dimension(fGhostNodes.Length(), the_field->NumDOF());
		fGhostNodesEquations.SetValueToPosition();
		fGhostNodesEquations += 1;
		
		fGhostIdToIndex.SetMap(fGhostNodes);
		fGhostIdToIndex.SetOutOfRange(InverseMapT::MinusOne);
	}

	/* clear values */
	K_G_NG = 0.0;

	/* form matrix */
	Particle().FormStiffness(fGhostIdToIndex, fGhostNodesEquations, K_G_NG);
}

/* set the field at the ghost nodes */
void FEManagerT_bridging::SetFieldValues(const StringT& field, const iArrayT& nodes, 
	const dArray2DT& values)
{
	const char caller[] = "FEManagerT_bridging::SetFieldValues";

#if __option(extended_errorcheck)
	if (nodes.Length() != values.MajorDim())
		ExceptionT::SizeMismatch(caller);
#endif

	/* get the associated field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* assume we're writing into the displacement array */
	dArray2DT& displacement = (*the_field)[0];
	for (int i = 0; i < values.MajorDim(); i++)	
		displacement.SetRow(nodes[i], values(i));

	/* reset the current configuration */
	fNodeManager->UpdateCurrentCoordinates();

	//NOTE: write the values into the KBC controller as well?
}

/* initialize nodes that follow the field computed by this instance */
void FEManagerT_bridging::InitInterpolation(const iArrayT& nodes, const StringT& field, 
	NodeManagerT& node_manager)
{
#pragma unused(field)

	fMainOut << "\n Number of interpolation points. . . . . . . . . = " << nodes.Length() << '\n';

	/* map nodes into cells (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().MaptoCells(nodes, &init_coords, NULL, fFollowerCellData);

	/* compute interpolation data */
	BridgingScale().InitInterpolation(nodes, fFollowerCellData);
}

/* field interpolations */
void FEManagerT_bridging::InterpolateField(const StringT& field, dArray2DT& nodal_values)
{
	/* interpolate in bridging scale element */
	BridgingScale().InterpolateField(field, fFollowerCellData, nodal_values);
}

/* return the interpolation matrix associated with the active degrees
 * of freedom */
void FEManagerT_bridging::InterpolationMatrix(const StringT& field, dSPMatrixT& G_Interpolation) const
{
	const char caller[] = "FEManagerT_bridging::InterpolationMatrix";

	/* get the associated field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* the shape functions values at the interpolating point */
	const dArray2DT& weights = fFollowerCellData.InterpolationWeights(); 

	/* redimension matrix if needed */
	int   ndof = the_field->NumDOF();
	int row_eq = weights.MajorDim()*ndof;	
	int col_eq = the_field->NumEquations();	
	if (G_Interpolation.Rows() != row_eq || G_Interpolation.Cols() != col_eq)
		G_Interpolation.Dimension(row_eq, col_eq, 0);

	/* clear */
	G_Interpolation = 0.0;

	/* element group information */
	const ContinuumElementT* continuum = fFollowerCellData.ContinuumElement();
	const iArrayT& cell = fFollowerCellData.InterpolatingCell();

	/* fill by rows - active DOF's only */
	row_eq = 0;
	for (int i = 0; i < weights.MajorDim(); i++)
	{
		/* element info */
		const ElementCardT& element_card = continuum->ElementCard(cell[i]);
		const iArrayT& eqnos = element_card.Equations();
		const iArrayT& nodes = element_card.NodesU();

		/* shape functions at the interpolation point */
		double* Na = weights(i);
		
		/* expand row dof's */
		for (int j = 0; j < ndof; j++) {
		
			/* expand col dof's */
			int eq_dex = 0;
			for (int k = 0; k < ndof; k++)
				for (int l = 0; l < nodes.Length(); l++) /* element nodes */
				{
					/* active value */
					int col_eq = eqnos[eq_dex++] - 1;
					if (col_eq > 0) /* write in */
						G_Interpolation.SetElement(row_eq, col_eq, Na[l]);
				}
		
			/* next row dof */
			row_eq++;
		}
	}
}

/* compute global interpolation matrix for all nodes whose support intersects MD region */
void FEManagerT_bridging::Ntf(dSPMatrixT& ntf, const iArrayT& atoms, iArrayT& activefenodes) const
{
	/* obtain global node numbers of nodes whose support intersects MD, create inverse map */
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();	// list of active nodes
	activefenodes = cell_nodes;
	const InverseMapT gtlnodes, gtlatoms;
	gtlnodes.SetMap(cell_nodes);	// create global to local map for active nodes
	gtlatoms.SetMap(atoms);		// create global to local map for all atoms - redundant?
	int numactivenodes = cell_nodes.Length();	// number of projected nodes
	int numatoms = atoms.Length();	// total number of atoms

	/* the shape functions values at the interpolating point */
	const dArray2DT& weights = fDrivenCellData.InterpolationWeights(); 

	/* redimension matrix if needed */
	int row_eq = numactivenodes;	// the number of projected nodes
	int col_eq = numatoms;	// total number of atoms
	if (ntf.Rows() != row_eq || ntf.Cols() != col_eq)
		ntf.Dimension(row_eq, col_eq, 0);

	/* clear */
	ntf = 0.0;
	
	/* element group information */
	const ContinuumElementT* continuum = fDrivenCellData.ContinuumElement();
	const iArrayT& cell = fDrivenCellData.InterpolatingCell();
	
	/* first loop over all atoms */
	for (int i = 0; i < col_eq; i++)
	{
		/* element info */
		const ElementCardT& element_card = continuum->ElementCard(cell[i]);
		const iArrayT& fenodes = element_card.NodesU();
		int dex = gtlatoms.Map(atoms[i]);	// global to local map for atoms
		
		/* put shape functions for nodes evaluated at each atom into global interpolation matrix */
		for (int j = 0; j < weights.MinorDim(); j++)
		{
			int dex2 = gtlnodes.Map(fenodes[j]);	// global to local map for nodes
			ntf.SetElement(dex2, dex, weights(dex,j));  // dex = i...
		}
	}
}

/* initialize data for the driving field */
void FEManagerT_bridging::InitProjection(const iArrayT& nodes, const StringT& field, 
	NodeManagerT& node_manager, bool make_inactive)
{
	const char caller[] = "FEManagerT_bridging::SetExactSolution";

	fMainOut << "\n Number of projection points . . . . . . . . . . = " << nodes.Length() << '\n';

	/* map nodes into cells (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().MaptoCells(nodes, &init_coords, NULL, fDrivenCellData);

	/* compute interpolation data */
	BridgingScale().InitInterpolation(nodes, fDrivenCellData);

	/* collect nodes in non-empty cells and generate cell connectivities 
	 * in local numbering*/
	fDrivenCellData.GenerateCellConnectivities();

	/* compute the mass matrix used for the projection */
	BridgingScale().InitProjection(fDrivenCellData);

	/* get the associated field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* create controller to driver solution */
	if (!fSolutionDriver) {
	
		/* construct new contoller */
		fSolutionDriver = new KBC_PrescribedT(*fNodeManager);
	
		/* add to field */
		the_field->AddKBCController(fSolutionDriver);
	}

	/* generate KBC cards - all degrees of freedom */
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();
	int ndof = the_field->NumDOF();
	if (make_inactive)
	{
		ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
		KBC_cards.Dimension(cell_nodes.Length()*ndof);
		int dex = 0;
		for (int j = 0; j < ndof; j++)
			for (int i = 0; i < cell_nodes.Length(); i++)
				KBC_cards[dex++].SetValues(cell_nodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	}

	/* dimension work space */
	fProjection.Dimension(cell_nodes.Length(), ndof);
	
	/* reset the group equations numbers */
	SetEquationSystem(the_field->Group());
}

/* project the point values onto the mesh */
void FEManagerT_bridging::ProjectField(const StringT& field, NodeManagerT& node_manager)
{
	const char caller[] = "FEManagerT_bridging::ProjectField";

	/* get the source field */
	FieldT* source_field = node_manager.Field(field);
	if (!source_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());

	/* compute the projection onto the mesh */
	const dArray2DT& source_field_values = (*source_field)[0];
	BridgingScale().ProjectField(field, fDrivenCellData, source_field_values, fProjection);

	/* write values into the field */
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();
	SetFieldValues(field, cell_nodes, fProjection);
}

/* project the point values onto the mesh */
void FEManagerT_bridging::InitialProject(const StringT& field, NodeManagerT& node_manager, dArray2DT& projectedu)
{
	const char caller[] = "FEManagerT_bridging::ProjectField";

	/* get the source field */
	FieldT* source_field = node_manager.Field(field);
	if (!source_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());

	/* compute the projection onto the mesh */
	const dArray2DT& source_field_values = (*source_field)[0];
	BridgingScale().InitialProject(field, fDrivenCellData, source_field_values, fProjection, projectedu);

	/* write values into the field */
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();
	SetFieldValues(field, cell_nodes, fProjection);
}

/* calculate the fine scale part of MD solution as well as total displacement u */
void FEManagerT_bridging::BridgingFields(const StringT& field, NodeManagerT& atom_node_manager, 
	NodeManagerT& fem_node_manager, dArray2DT& totalu)
{
	const char caller[] = "FEManagerT_bridging::ProjectField";

	/* get the fem and md fields */
	FieldT* atom_field = atom_node_manager.Field(field);
	FieldT* fem_field = fem_node_manager.Field(field);
	if (!atom_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());
	if (!fem_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());
	
	/* compute the fine scale part of MD solution as well as total displacement u */
	const dArray2DT& atom_values = (*atom_field)[0];
	const dArray2DT& fem_values = (*fem_field)[0];
	BridgingScale().BridgingFields(field, fDrivenCellData, atom_values, fem_values, fProjection, totalu);
}


/* the residual for the given group */
const dArrayT& FEManagerT_bridging::Residual(int group) const 
{
	return fSolvers[group]->RHS(); 
}

/* set the reference error for the given group */
void FEManagerT_bridging::SetReferenceError(int group, double error) const
{
	/* retrieve nonlinear solver */
	NLSolver* solver = dynamic_cast<NLSolver*>(fSolvers[group]);

	/* silent in failuer */
	if (solver) solver->SetReferenceError(error);
}

/* return the internal forces for the given solver group associated with the
 * most recent call to FEManagerT_bridging::FormRHS. */
const dArray2DT& FEManagerT_bridging::InternalForce(int group) const
{
//TEMP - only look for the contribution from the particle group
	const char caller[] = "FEManagerT_bridging::InternalForce";
	if (!fParticle)
		ExceptionT::GeneralFail(caller, "particle group not set");
	else if (!fParticle->InGroup(group))
		ExceptionT::GeneralFail(caller, "particles are not in group %d", group);
		
	return fParticle->InternalForce(group);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* initialize solver information */
void FEManagerT_bridging::SetSolver(void)
{
	/* dimension list of cumulative update vectors */
	fCumulativeUpdate.Dimension(NumGroups());
	
	/* dimension list of pointers to external force vectors */
	fExternalForce.Dimension(NumGroups());
	fExternalForce = NULL;

	fExternalForce2D.Dimension(NumGroups());
	fExternalForce2DNodes.Dimension(NumGroups());
	fExternalForce2DEquations.Dimension(NumGroups());
	fExternalForce2D = NULL;
	fExternalForce2DNodes = NULL;

	/* inherited */
	FEManagerT::SetSolver();
}

/*************************************************************************
 * Private
 *************************************************************************/

/* the particle group */
ParticleT& FEManagerT_bridging::Particle(void) const
{
	/* find bridging scale group */
	if (!fParticle) {
	
		/* search through element groups */
		for (int i = 0; !fParticle && i < fElementGroups.Length(); i++)
		{
			/* try cast */
			ElementBaseT* element_base = fElementGroups[i];
			
			/* need non-const pointer to this */
			FEManagerT_bridging* fe = (FEManagerT_bridging*) this;
			fe->fParticle = dynamic_cast<ParticleT*>(element_base);
		}
		
		/* not found */
		if (!fParticle)
			ExceptionT::GeneralFail("FEManagerT_bridging::Particle",
				"did not find ParticleT group");
	}
	
	return *fParticle;
}

/* the bridging scale element group */
BridgingScaleT& FEManagerT_bridging::BridgingScale(void) const
{
	/* find bridging scale group */
	if (!fBridgingScale) {
	
		/* search through element groups */
		for (int i = 0; !fBridgingScale && i < fElementGroups.Length(); i++)
		{
			/* try cast */
			ElementBaseT* element_base = fElementGroups[i];
			
			/* need non-const pointer to this */
			FEManagerT_bridging* fe = (FEManagerT_bridging*) this;
			fe->fBridgingScale = dynamic_cast<BridgingScaleT*>(element_base);
		}
		
		/* not found */
		if (!fBridgingScale)
			ExceptionT::GeneralFail("FEManagerT_bridging::BridgingScale",
				"did not find BridgingScaleT element group");
	}
	
	return *fBridgingScale;
}

#endif  /* BRIDGING_ELEMENT */
