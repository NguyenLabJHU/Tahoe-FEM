/* $Id: FEManagerT_bridging.cpp,v 1.16.4.3 2004-04-03 03:18:39 paklein Exp $ */
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
#include "ParticlePairT.h"
#include "dSPMatrixT.h"

/* headers needed to compute the correction for overlap */
#include "SolidMatListT.h"
#include "FCC3D.h"
#include "Hex2D.h"
#include "BondLatticeT.h"

using namespace Tahoe;

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT(input, output, comm),
	fBridgingIn(bridging_input),
	fBridgingScale(NULL),
	fSolutionDriver(NULL)
{
/* class requires RTTI */
#ifdef __NO_RTTI__
#pragma message("requires RTTI")
ExceptionT::GeneralFail("FEManagerT_bridging::FEManagerT_bridging", "requires RTTI");
#endif
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

/* compute internal correction for the overlap region */
void FEManagerT_bridging::CorrectOverlap(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap";

	/* collect nodes and cells in the overlap region */
	iArrayT overlap_cell;
	iArrayT overlap_node;
	CollectOverlapRegion(overlap_cell, overlap_node);
	if (overlap_node.Length() == 0) return;

	/* Cauchy-Born constitutive model */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	const MaterialListT& mat_list = coarse->MaterialsList();
	if (mat_list.Length() > 2) ExceptionT::GeneralFail(caller, "expecting at most 2 materials, not %d", mat_list.Length());

	ContinuumMaterialT* cont_mat = mat_list[0];
	Hex2D* hex_2D = dynamic_cast<Hex2D*>(cont_mat);
	FCC3D* fcc_3D = dynamic_cast<FCC3D*>(cont_mat);
	if (!hex_2D && !fcc_3D) ExceptionT::GeneralFail(caller, "could not resolve C-B material");

	const BondLatticeT& bond_lattice = (hex_2D) ? hex_2D->BondLattice() : fcc_3D->BondLattice();
	const dArray2DT& bonds = bond_lattice.Bonds();
	double cell_volume = (hex_2D) ? hex_2D->CellVolume() : fcc_3D->CellVolume();

	/* initialize unknowns */
	dArray2DT bond_densities(overlap_cell.Length(), coarse->NumIP()*bonds.MajorDim());
	bond_densities = 1.0;

	/* solve unkowns */

	/* change overlap element material numbers */

	/* write unknowns into the state variable space */
	
}

/* (re-)set the equation number for the given group */
void FEManagerT_bridging::SetEquationSystem(int group, int start_eq_shift)
{
	/* inherited */
	FEManagerT::SetEquationSystem(group, start_eq_shift);

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
		ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
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
	ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
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

/* compute the product with transpose of the interpolation matrix */
void FEManagerT_bridging::MultNTf(const PointInCellDataT& N, const dArray2DT& f, const iArrayT& f_rows, dArray2DT& NTf) const
{
	/* coarse scale element group */
	const ContinuumElementT* continuum = N.ContinuumElement();

	/* interpolation "matrix" */
	const InverseMapT& f_rows_map = N.GlobalToLocal();
	const dArray2DT& interpolation_weights = N.InterpolationWeights(); 
	const iArrayT& cell = N.InterpolatingCell();
	
	/* loop over active rows in f and accumulate */
	dArrayT weights;
	dArrayT f_row;
	dArrayT NTf_row;
	for (int i = 0; i < f_rows.Length(); i++) {
	
		/* mappings */
		int f_row_i = f_rows[i];
		int f_row_i_map = f_rows_map.Map(f_row_i);
		f.RowAlias(f_row_i, f_row);
		
		/* cell nodes */
		const iArrayT& NTf_rows = continuum->ElementCard(cell[f_row_i_map]).NodesU();
		
		/* interpolation weights */
		interpolation_weights.RowAlias(f_row_i_map, weights);
		
		/* loop over neighbors */
		for (int j = 0; j < weights.Length(); j++) {
		
			/* row in output */
			NTf.RowAlias(NTf_rows[j], NTf_row);
		
			/* accumulate contribution */
			NTf_row.AddScaled(weights[j], f_row);
		}
	}
}

void FEManagerT_bridging::MultNTf(const InterpolationDataT& N, const dArray2DT& f, const iArrayT& f_rows, dArray2DT& NTf) const
{
	/* interpolation "matrix" */
	const InverseMapT& f_rows_map = N.Map();
	const RaggedArray2DT<double>& neighbor_weights = N.NeighborWeights(); 
	const RaggedArray2DT<int>& neighbors = N.Neighbors();
	
	/* loop over active rows in f and accumulate */
	dArrayT weights;
	dArrayT f_row;
	iArrayT NTf_rows;
	dArrayT NTf_row;
	for (int i = 0; i < f_rows.Length(); i++) {
	
		/* mappings */
		int f_row_i = f_rows[i];
		int f_row_i_map = f_rows_map.Map(f_row_i);
		f.RowAlias(f_row_i, f_row);
		
		/* neigbors */
		neighbors.RowAlias(f_row_i_map, NTf_rows);
		
		/* interpolation weights */
		neighbor_weights.RowAlias(f_row_i_map, weights);
		
		/* loop over neighbors */
		for (int j = 0; j < weights.Length(); j++) {
		
			/* row in output */
			NTf.RowAlias(NTf_rows[j], NTf_row);
		
			/* accumulate contribution */
			NTf_row.AddScaled(weights[j], f_row);
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

/* transpose follower cell data */
void FEManagerT_bridging::TransposeFollowerCellData(InterpolationDataT& transpose)
{
	/* map from global point id to row in interpolation data */
	const InverseMapT& map = fFollowerCellData.GlobalToLocal();

	/* interpolation weights */
	const dArray2DT& neighbor_weights = fFollowerCellData.InterpolationWeights();

	/* collect connectivities */
	iArray2DT neighbors(neighbor_weights.MajorDim(), neighbor_weights.MinorDim());
	const iArrayT& cell = fFollowerCellData.InterpolatingCell();
	const ContinuumElementT* continuum = fFollowerCellData.ContinuumElement();
	for (int i = 0; i < neighbors.MajorDim(); i++) {

		/* element nodes */
		int element = cell[i];
		const iArrayT& nodes = continuum->ElementCard(element).NodesU();

		/* copy */
		neighbors.SetRow(i, nodes);
	}

	transpose.Transpose(map, neighbors, neighbor_weights);
}

/* set the reference error for the given group */
void FEManagerT_bridging::SetReferenceError(int group, double error) const
{
	/* retrieve nonlinear solver */
	NLSolver* solver = dynamic_cast<NLSolver*>(fSolvers[group]);

	/* silent in failure */
	if (solver) solver->SetReferenceError(error);
}

/* return the internal forces for the given solver group associated with the
 * most recent call to FEManagerT_bridging::FormRHS. */
const dArray2DT& FEManagerT_bridging::InternalForce(int group) const
{
	const char caller[] = "FEManagerT_bridging::InternalForce";

	/* search through element groups */
	ElementBaseT* element = NULL;
	for (int i = 0; i < fElementGroups->Length(); i++) {
		ElementBaseT* element_tmp = (*fElementGroups)[i];
		if (element_tmp != fBridgingScale && element_tmp->InGroup(group))
		{
			/* already found element group */
			if (element) ExceptionT::GeneralFail(caller, "solver group %d contains more than one element group", group);
			element = element_tmp;
		}
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
	ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
	if (!particle)
		ExceptionT::GeneralFail("FEManagerT_bridging::PropertiesMap",
			"group %d is not a particle group", element_group);
	return particle->PropertiesMap();
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

void FEManagerT_bridging::CollectOverlapRegion(iArrayT& overlap_cell, iArrayT& overlap_node) const
{
	const char caller[] = "FEManagerT_bridging::CollectOverlapRegion";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");

	/* dimensions */
	int nnd = fNodeManager->NumNodes();
	int nel = coarse->NumElements();
	int nen = coarse->NumElementNodes();

	/* mark nodes that aren't active */
	ArrayT<char> is_overlap_node(nnd);
	is_overlap_node = 't';
	const iArrayT& projected_nodes = fDrivenCellData.CellNodes();
	for (int i = 0; i < projected_nodes.Length(); i++) 
		is_overlap_node[projected_nodes[i]] = 'f';

	/* find cells in overlap region */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();
	ArrayT<char> is_overlap_cell(nel);
	is_overlap_cell = 'f';
	int num_overlap_cell = 0;
	iArrayT cell_point;
	for (int i = 0; i < nel; i++) /* cells must contain follower points bonded to free points */
		if (point_in_cell.MinorDim(i) > 0) {
			point_in_cell.RowAlias(i, cell_point);
			bool has_follower = false;
			for (int j = 0; !has_follower && j < cell_point.Length(); j++) /* find a follower point */
				if (follower_point_map.Map(cell_point[j]) > -1) {
					has_follower = true;
					is_overlap_cell[i] = 't';
					num_overlap_cell++;
				}		
		}
	
	for (int i = 0; i < nel; i++) /* cells must contain at least one active node */
		if (is_overlap_cell[i] == 't') {
		
			/* check cell nodes */
			const iArrayT& nodes = coarse->ElementCard(i).NodesU();
			bool OK = false;
			for (int j = 0; !OK && j < nen; j++)
				if (is_overlap_node[nodes[j]] == 't')
					OK = true;
			
			/* remove cell */
			if (!OK) {
				is_overlap_cell[i] = 'f';
				num_overlap_cell--;			
			}	
		}
		
	overlap_cell.Dimension(num_overlap_cell);
	num_overlap_cell = 0;
	for (int i = 0; i < nel; i++)
		if (is_overlap_cell[i] == 't')
			overlap_cell[num_overlap_cell++] = i;
			
	if (fPrintInput) {
		overlap_cell++;
		fMainOut << "\n overlap cells: " << overlap_cell.Length() << '\n';
		fMainOut << overlap_cell.wrap(5) << '\n';
		overlap_cell--;	
	}
	
	/* find nodes in overlap region */
	is_overlap_node = 'f';
	num_overlap_cell = overlap_cell.Length();
	int num_overlap_node = 0;
	for (int i = 0; i < num_overlap_cell; i++) /* support must include overlap cells */ {
		const iArrayT& nodes = coarse->ElementCard(overlap_cell[i]).NodesU();
		for (int j = 0; j < nen; j++) {
			char& t_f = is_overlap_node[nodes[j]];
			if (t_f == 'f') {
				t_f = 't';
				num_overlap_node++;
			}		
		}	
	}
	
	for (int i = 0; i < projected_nodes.Length(); i++) /* must be a free node */ {
		char& t_f = is_overlap_node[projected_nodes[i]];
		if (t_f == 't') /* remove node */ {
			t_f = 'f';
			num_overlap_node--;
		}
	}

	overlap_node.Dimension(num_overlap_node);
	num_overlap_node = 0;
	for (int i = 0; i < nnd; i++)
		if (is_overlap_node[i] == 't')
			overlap_node[num_overlap_node++] = i;

	if (fPrintInput) {
		overlap_node++;
		fMainOut << "\n overlap nodes: " << overlap_node.Length() << '\n';
		fMainOut << overlap_node.wrap(5) << '\n';
		overlap_node--;	
	}	
}

/* return the given instance of the ParticlePairT element group or NULL if not found */
const ParticlePairT* FEManagerT_bridging::ParticlePair(int instance) const
{
	/* search through element groups for particles */
	int count = 0;
	for (int i = 0; i < fElementGroups->Length(); i++)
	{
		/* pointer to element group */
		ElementBaseT* element_base = (*fElementGroups)[i];
		
		/* attempt cast to particle type */
		ParticlePairT* particle_pair = dynamic_cast<ParticlePairT*>(particle_pair);
		if (particle_pair && count++ == instance)
			return particle_pair;
	}

	/* fall through */
	return NULL;
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
