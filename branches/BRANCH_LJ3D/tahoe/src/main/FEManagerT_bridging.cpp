/* $Id: FEManagerT_bridging.cpp,v 1.19.2.1 2004-04-10 20:56:38 hspark Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "KBC_PrescribedT.h"
#include "KBC_CardT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "NLSolver.h"
#include "CommManagerT.h"

#include "BridgingScaleT.h"
#include "ParticleT.h"
#include "dSPMatrixT.h"
#include "EAMFCC3D.h"
#include "EAMT.h"
#include "ShapeFunctionT.h"
#include "dSymMatrixT.h"
#include "ElementSupportT.h"

using namespace Tahoe;

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT(input, output, comm),
	fBridgingIn(bridging_input),
	fBridgingScale(NULL),
	fSolutionDriver(NULL),
	fEAMFCC3D(NULL),
	fEAMT(NULL)
{

}

/* destructor */
FEManagerT_bridging::~FEManagerT_bridging(void)
{
	delete fEAMFCC3D;
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
	eqnos.Dimension(activefenodes.Length(), thefield->NumDOF());
	thefield->SetLocalEqnos(activefenodes, eqnos);	// crashes here if not nodes and atoms everywhere
}

/* initialize the ghost node information */
void FEManagerT_bridging::InitGhostNodes(bool include_image_nodes)
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

	/* search through element groups for particles */
	bool found = false;
	for (int i = 0; i < fElementGroups->Length(); i++)
	{
		/* pointer to element group */
		ElementBaseT* element_base = (*fElementGroups)[i];
		
		/* attempt cast to particle type */
#ifndef __NO_RTTI_
		ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
#else /* no RTTI */
		ParticleT* particle = element_base->dynamic_cast_ParticleT();
#endif
		if (particle) 
		{
			found = true;
			particle->SetSkipParticles(fGhostNodes);
			particle->SetConfiguration();
		}
	}
	if (!found) ExceptionT::GeneralFail(caller, "no particle group found");
	
	/* reset the group equations numbers */
	SetEquationSystem(the_field->Group());

	/* echo ghost nodes */
	if (fPrintInput) {
		fGhostNodes++;
		fMainOut << "\n Ghost nodes:\n";
		fMainOut << fGhostNodes.wrap(5) << '\n';
		fGhostNodes--;
	}

	/* initialize potential non-ghost nodes */
	CommManagerT* comm = FEManagerT::CommManager();	
	const ArrayT<int>* part_nodes = comm->PartitionNodes();
	iArrayT is_ghost;
	if (include_image_nodes || !part_nodes) {
		/* assuming there are no images in the list of ghost nodes */
		fNonGhostNodes.Dimension(fModelManager->NumNodes() - fGhostNodes.Length());
		is_ghost.Dimension(fModelManager->NumNodes());
		is_ghost = 0;	
	} else { /* remove image nodes */
		is_ghost.Dimension(fModelManager->NumNodes());
		is_ghost = 1;

		/* initialize potential non-ghost nodes */		
		const int* p = part_nodes->Pointer();
		int npn = part_nodes->Length();
		for (int i = 0; i < npn; i++)
			is_ghost[*p++] = 0;
	}	

	/* mark nodes as ghost */
	for (int i = 0; i < fGhostNodes.Length(); i++) {
		int& is_ghost_i = is_ghost[fGhostNodes[i]];
		if (is_ghost_i == 1)
			ExceptionT::GeneralFail(caller, "ghost node %d is duplicated or image",
				fGhostNodes[i]+1);
		else
			is_ghost_i = 1;
	}

	/* collect non-ghost nodes */
	if (fNonGhostNodes.Length() == 0) 
		fNonGhostNodes.Dimension(is_ghost.Count(0));
	dex = 0;
	for (int i = 0; i < is_ghost.Length(); i++)
		if (is_ghost[i] == 0)
			fNonGhostNodes[dex++] = i;
}

/* prescribe the motion of ghost nodes */
void FEManagerT_bridging::SetGhostNodeKBC(KBC_CardT::CodeT code, const dArray2DT& values)
{
	const char caller[] = "FEManagerT_bridging::SetGhostNodeKBC";
	if (!fSolutionDriver) ExceptionT::GeneralFail(caller, "controller for ghost node motion not set");

	/* fetch cards */
	ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();

	/* check dimensions */
	int ndof = values.MinorDim();
	if (KBC_cards.Length()/ndof != values.MajorDim())
		ExceptionT::SizeMismatch(caller, "expecting %d nodal values not %d",
			KBC_cards.Length()/ndof, values.MajorDim());

	/* loop over cards */
	for (int i = 0; i < KBC_cards.Length(); i++)
	{
		/* retrieve values set during InitGhostNodes */
		KBC_CardT& card = KBC_cards[i];
		int node = card.Node();
		int dof  = card.DOF();
		int schd = card.ScheduleNum();
	
		/* reset code and value */
		card.SetValues(node, dof, code, schd, values[i]);
	}
}

/* compute the ghost-nonghost part of the stiffness matrix */
void FEManagerT_bridging::Form_G_NG_Stiffness(const StringT& field, int element_group, dSPMatrixT& K_G_NG)
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

	/* try cast */
	ElementBaseT* element_base = (*fElementGroups)[element_group];
#ifndef __NO_RTTI_
	ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
#else
	ParticleT* particle = element_base->dynamic_cast_ParticleT();
#endif
	if (!particle) ExceptionT::GeneralFail(caller, "element group %d is not a particle group", element_group);

	/* form matrix */
	particle->FormStiffness(fGhostIdToIndex, fGhostNodesEquations, K_G_NG);
}

/* set the field at the ghost nodes */
void FEManagerT_bridging::SetFieldValues(const StringT& field, const iArrayT& nodes, int order, 
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

	/* can write into any field due to order parameter */
	dArray2DT& arbitraryfield = (*the_field)[order];
	for (int i = 0; i < values.MajorDim(); i++)	
		arbitraryfield.SetRow(nodes[i], values(i));

	/* reset the current configuration */
	fNodeManager->UpdateCurrentCoordinates();

	//NOTE: write the values into the KBC controller as well?
}

/* return the "lumped" (scalar) mass associated with the given nodes */
void FEManagerT_bridging::LumpedMass(const iArrayT& nodes, dArrayT& mass) const
{
	/* initialize */
	mass.Dimension(nodes.Length());
	mass = 0.0;

	/* accumulate element contribution */
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->LumpedMass(nodes, mass);
}

/* initialize nodes that follow the field computed by this instance */
void FEManagerT_bridging::InitInterpolation(const iArrayT& nodes, const StringT& field, 
	NodeManagerT& node_manager)
{
#pragma unused(field)

	fMainOut << "\n Number of interpolation points. . . . . . . . . = " << nodes.Length() << '\n';

	/* compute interpolation data (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().InitInterpolation(nodes, &init_coords, NULL, fFollowerCellData);
}

/* field interpolations */
void FEManagerT_bridging::InterpolateField(const StringT& field, int order, dArray2DT& nodal_values)
{
	/* interpolate in bridging scale element */
	BridgingScale().InterpolateField(field, order, fFollowerCellData, nodal_values);
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
		const double* Na = weights(i);
		
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
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();	// list of active nodes fDrivenCellData
	activefenodes = cell_nodes;
	InverseMapT gtlnodes;
	gtlnodes.SetMap(cell_nodes);	// create global to local map for active nodes
	int numactivenodes = cell_nodes.Length();	// number of projected nodes
	int numatoms = atoms.Length();	// total number of non ghost atoms

	/* the shape functions values at the interpolating point */
	const dArray2DT& weights = fDrivenCellData.InterpolationWeights(); 

	/* dimension matrix if needed */
	int row_eq = numactivenodes;	// the number of projected nodes
	int col_eq = numatoms;	// total number of non ghost atoms
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

		/* put shape functions for nodes evaluated at each atom into global interpolation matrix */
		for (int j = 0; j < weights.MinorDim(); j++)
		{
			int dex2 = gtlnodes.Map(fenodes[j]);	// global to local map for nodes
			ntf.SetElement(dex2, i, weights(i,j));  // dex = i...
		}
	}
}

/* initialize data for the driving field */
void FEManagerT_bridging::InitProjection(CommManagerT& comm, const iArrayT& nodes, const StringT& field, 
	NodeManagerT& node_manager, bool make_inactive)
{
	const char caller[] = "FEManagerT_bridging::InitProjection";
	fMainOut << "\n Number of projection points . . . . . . . . . . = " << nodes.Length() << '\n';

	/* initialize the projection (using reference coordinates) */
	const dArray2DT& init_coords = node_manager.InitialCoordinates();
	BridgingScale().InitProjection(comm, nodes, &init_coords, NULL, fDrivenCellData);

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

	/* construct EAMFCC3D pointer if 3D bridging scale using EAM */
	//if (the_field->NumDOF() == 3)
	//{
		/* construct EAMFCC3D pointer if 3D bridging scale using EAM */
	//	ifstreamT& in = Input();
	//	fEAMFCC3D = new EAMFCC3D(in, 4, 3, 54);
	//	fEAMFCC3D->InitBondTables();
	//}
}

/* indicate whether image nodes should be included in the projection */
bool FEManagerT_bridging::ProjectImagePoints(void) const
{
	return BridgingScale().ProjectImagePoints();
}

/* project the point values onto the mesh */
void FEManagerT_bridging::ProjectField(const StringT& field, const NodeManagerT& node_manager, int order)
{
	const char caller[] = "FEManagerT_bridging::ProjectField";

	/* get the source field */
	const FieldT* source_field = node_manager.Field(field);
	if (!source_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());

	/* compute the projection onto the mesh */
	const dArray2DT& source_field_values = (*source_field)[0];
	BridgingScale().ProjectField(fDrivenCellData, source_field_values, fProjection);

	/* write values into the field */
	const iArrayT& cell_nodes = fDrivenCellData.CellNodes();
	SetFieldValues(field, cell_nodes, order, fProjection);
}

/* compute the coarse scale projection at the source points */
void FEManagerT_bridging::CoarseField(const StringT& field, const NodeManagerT& node_manager, int order, 
	dArray2DT& coarse)
{
	const char caller[] = "FEManagerT_bridging::ProjectField";

	/* get the source field */
	const FieldT* source_field = node_manager.Field(field);
	if (!source_field) ExceptionT::GeneralFail(caller, "could not resolve source field \"%s\"", field.Pointer());
	const dArray2DT& field_values = (*source_field)[order];

	/** compute the coarse scale part of the source field */
	BridgingScale().CoarseField(fDrivenCellData, field_values, coarse);
}

/* project the point values onto the mesh */
void FEManagerT_bridging::InitialProject(const StringT& field, NodeManagerT& node_manager, dArray2DT& projectedu,
int order)
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
	SetFieldValues(field, cell_nodes, order, fProjection);
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

/* set the reference error for the given group */
void FEManagerT_bridging::SetReferenceError(int group, double error) const
{
	/* retrieve nonlinear solver */
  NLSolver* solver = TB_DYNAMIC_CAST(NLSolver*, fSolvers[group]);

	/* silent in failuer */
	if (solver) solver->SetReferenceError(error);
}

/* return the internal forces for the given solver group associated with the
 * most recent call to FEManagerT_bridging::FormRHS. */
const dArray2DT& FEManagerT_bridging::InternalForce(int group) const
{
	const char caller[] = "FEManagerT_bridging::InternalForce";

	/* search through element groups */
	ElementBaseT* element = NULL;
	for (int i = 0; i < fElementGroups->Length(); i++)
		if ((*fElementGroups)[i]->InGroup(group))
		{
			/* already found element group */
			if (element) ExceptionT::GeneralFail(caller, "solver group %d contains more than one element group", group);
			element = (*fElementGroups)[i];
		}
	
	/* no elements in the group */
	if (!element) ExceptionT::GeneralFail(caller, "no elements in solver group %d", group);

	return element->InternalForce(group);
}

/* return the properties map for the given element group */
nMatrixT<int>& FEManagerT_bridging::PropertiesMap(int element_group)
{
	/* try cast to particle type */
	ElementBaseT* element_base = (*fElementGroups)[element_group];
#ifndef __NO_RTTI__
	ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
#else
	ParticleT* particle = element_base->dynamic_cast_ParticleT();
#endif
	if (!particle)
		ExceptionT::GeneralFail("FEManagerT_bridging::PropertiesMap",
			"group %d is not a particle group", element_group);
	return particle->PropertiesMap();
}

/* calculate EAM total electron density at ghost atoms */
void FEManagerT_bridging::ElecDensity(int length, dArray2DT& elecdens, dArray2DT& embforce)
{
	/* try constructing an EAMFCC3D here - need to construct new EAMFCC3D for each ghost atom? */
	StringT field = "displacement";
	const ContinuumElementT* continuum = fFollowerCellData.ContinuumElement();
	const ElementCardT& element_card1 = continuum->ElementCard(0);
	const iArrayT& nodes1 = element_card1.NodesU();
	int nen = nodes1.Length();	// number of nodes per element
	const ShapeFunctionT& shape = continuum->ShapeFunction();	// access element shape functions

	/* get the field */
	const FieldT* the_field = fNodeManager->Field(field);
	LocalArrayT loc_field(LocalArrayT::kDisp, nen, the_field->NumDOF());
	int nsd = the_field->NumDOF();
	loc_field.SetGlobal((*the_field)[0]); /* displacement only */
	
	const iArrayT& cell = fFollowerCellData.InterpolatingCell();
	RaggedArray2DT<double>& inversemap = fFollowerCellData.PointInCellCoords(); 
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& global_to_local = fFollowerCellData.GlobalToLocal();
	dArrayT Na, sourcea, sourceb(nsd);
	dMatrixT fgrad(nsd), eye(nsd), green(nsd);
	eye.Identity(1.0);
	dArray2DT DNa;
	dSymMatrixT green1(dSymMatrixT::k3D);
	double ed, ef;
	
	/* loop over all elements which contain atoms/ghost atoms */
	for (int i = 0; i < inversemap.MajorDim(); i++)
	{
		int np = inversemap.MinorDim(i);
		int np1 = point_in_cell.MinorDim(i);
		if (np > 0)
		{
			const int* points = point_in_cell(i);	// global # of atom in cell
			const double* inverse = inversemap(i);	// parent domain shape function values
			inversemap.RowAlias(i, sourcea);
			for (int j = 0; j < np1; j++)
			{
				int point_dex = global_to_local.Map(points[j]);
				
				/* because ghost atoms stored first in boundaryghost atoms /*
				/* global to local map has ghost atoms first, boundary atoms second */
				/* this loop ignores the boundary atoms */
				if (point_dex < length)
				{
					/* determine the element the ghost atom lies in */
					const ElementCardT& element_card = continuum->ElementCard(cell[point_dex]);
					const iArrayT& nodes = element_card.NodesU();

					/* collect local displacements */
					loc_field.SetLocal(nodes);
			
					/* obtain parent coordinates */
					sourceb.CopyPart(0, sourcea, j*nsd, nsd);
				
					/* calculate gradient of displacement field */
					shape.GradU(loc_field, fgrad, sourceb, Na, DNa);
				
					/* calculate deformation gradient = 1 + GradU */
					fgrad+=eye;
	
					/* calculate green strain = .5*(F^{T}F-I) */
					green.MultATB(fgrad, fgrad, 0);
					green-=eye;
					green*=.5;
					green1.Symmetrize(green);
				
					/* calculate/store electron density/embedding force for each ghost atom */
					fEAMFCC3D->ElectronDensity(green1, ed, ef);
					elecdens(point_dex, 0) = ed;
					embforce(point_dex, 0) = ef;
				}
			}
		}
	}
}

/* add external electron density contribution to ghost atoms */
void FEManagerT_bridging::SetExternalElecDensity(const dArray2DT& elecdens, const iArrayT& ghostatoms)
{
	const char caller[] = "FEManagerT_bridging::SetExternalElecDensity";

	/* check */
	if (ghostatoms.Length() != elecdens.MajorDim()) 
		ExceptionT::SizeMismatch(caller);

	/* store pointers in EAMT */
	EAM().SetExternalElecDensity(elecdens, ghostatoms);
}

/* add external embedding force contribution to ghost atoms */
void FEManagerT_bridging::SetExternalEmbedForce(const dArray2DT& embforce, const iArrayT& ghostatoms)
{
	const char caller[] = "FEManagerT_bridging::SetExternalForce";

	/* check */
	if (ghostatoms.Length() != embforce.MajorDim()) 
		ExceptionT::SizeMismatch(caller);

	/* store pointers in EAMT */
	EAM().SetExternalEmbedForce(embforce, ghostatoms);
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

/* the bridging scale element group */
BridgingScaleT& FEManagerT_bridging::BridgingScale(void) const
{
	/* find bridging scale group */
	if (!fBridgingScale) {
	
		/* search through element groups */
		for (int i = 0; !fBridgingScale && i < fElementGroups->Length(); i++)
		{
			/* try cast */
			ElementBaseT* element_base = (*fElementGroups)[i];
			
			/* need non-const pointer to this */
			FEManagerT_bridging* fe = (FEManagerT_bridging*) this;
#ifndef __NO_RTTI__
			fe->fBridgingScale = dynamic_cast<BridgingScaleT*>(element_base);
#else
			fe->fBridgingScale = element_base->dynamic_cast_BridgingScaleT();
#endif
		}
		
		/* not found */
		if (!fBridgingScale)
			ExceptionT::GeneralFail("FEManagerT_bridging::BridgingScale",
				"did not find BridgingScaleT element group");
	}
	
	return *fBridgingScale;
}

/* the EAMT element group */
EAMT& FEManagerT_bridging::EAM(void) const
{
	/* find EAMT group */
	if (!fEAMT) 
	{
	
		/* search through element groups */
		for (int i = 0; !fEAMT && i < fElementGroups->Length(); i++)
		{
			/* try cast */
			ElementBaseT* element_base = (*fElementGroups)[i];
			
			/* need non-const pointer to this */
			FEManagerT_bridging* fe = (FEManagerT_bridging*) this;
			fe->fEAMT = dynamic_cast<EAMT*>(element_base);
		}
		
		/* not found */
		if (!fEAMT)
			ExceptionT::GeneralFail("FEManagerT_bridging::EAM",
				"did not find EAMT element group");
	}
	
	return *fEAMT;
}

#endif  /* BRIDGING_ELEMENT */