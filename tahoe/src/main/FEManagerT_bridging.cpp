/* $Id: FEManagerT_bridging.cpp,v 1.34 2005-04-07 19:18:05 d-farrell2 Exp $ */

#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "KBC_ControllerT.h"
#include "KBC_CardT.h"
#include "NLSolver.h"
#include "CommManagerT.h"
#include "BridgingScaleT.h"
#include "ParticleT.h"
#include "ParticlePairT.h"
#include "dSPMatrixT.h"
#include "EAMFCC3D.h"
#include "EAMT.h"
#include "ShapeFunctionT.h"
#include "dSymMatrixT.h"
#include "ElementSupportT.h"

/* headers needed to compute the correction for overlap */
#include "ContinuumElementT.h"
#include "SolidMatListT.h"
#include "FCC3D.h"
#include "Hex2D.h"
#include "Chain1D.h"
#include "BondLatticeT.h"
#include "nArrayGroupT.h"
#include "nVariMatrixT.h"

#include "LAdMatrixT.h" //TEMP
#include "CCSMatrixT.h"

/* debugging */
#define __DEBUG__ 1

/* atom/point types */
const char free_ = 'f';
const char not_free_ = 'n';

/* element types */
const char p_0 = 'a'; /* bond density = 0 */
const char p_1 = 'b'; /* bond density = 1 */
const char p_x = 'c'; /* unknown: 0 < bond density < 1 */

using namespace Tahoe;

/* constructor */
FEManagerT_bridging::FEManagerT_bridging(const StringT& input, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	FEManagerT(input, output, comm, argv, task),
	fBridgingScale(NULL),
	fSolutionDriver(NULL),
	fEAMFCC3D(NULL),
	fEAMT(NULL)
{
	SetName("tahoe_bridging");
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

/* enforce zero bond density in projected cells */
void FEManagerT_bridging::DeactivateFollowerCells(void)
{
	const char caller[] = "FEManagerT_bridging::DeactivateFollowerCells";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fDrivenCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "projection data not set");
	const MaterialListT& mat_list = coarse->MaterialsList();
	if (mat_list.Length() > 1) 
		ExceptionT::GeneralFail(caller, "expecting only 1 material not %d", mat_list.Length());

	ContinuumMaterialT* cont_mat = mat_list[0];
	Hex2D* hex_2D = dynamic_cast<Hex2D*>(cont_mat);
	FCC3D* fcc_3D = dynamic_cast<FCC3D*>(cont_mat);
	Chain1D* mat_1D = dynamic_cast<Chain1D*>(cont_mat);
	if (!mat_1D && !hex_2D && !fcc_3D) ExceptionT::GeneralFail(caller, "could not resolve C-B material");

	const BondLatticeT& bond_lattice = (hex_2D) ? hex_2D->BondLattice() : ((fcc_3D) ? fcc_3D->BondLattice() : mat_1D->BondLattice());
	const dArray2DT& bonds = bond_lattice.Bonds();
	int num_densities = coarse->NumIP()*bonds.MajorDim();	

	/* collect cells in projected region */
	iArrayT cells;
	BridgingScale().CollectProjectedCells(fDrivenCellData, cells);

	/* write unknowns into the state variable space */
	ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
	for (int i = 0; i < cells.Length(); i++) {
	
		/* element information */
		ElementCardT& element = non_const_coarse->ElementCard(cells[i]);
	
		/* allocate space */
		element.Dimension(0, num_densities);
		
		/* zero bond densities */
		element.DoubleData() = 0.0;
	}
}

/* (re-)set the equation number for the given group */
void FEManagerT_bridging::SetEquationSystem(int group, int start_eq_shift)
{
	/* inherited */
	FEManagerT::SetEquationSystem(group, start_eq_shift);

	//NOTE: this is going to break if the equation numbers has changed since the force was set
//	if (fExternalForce2D[group])
//		ExceptionT::GeneralFail("FEManagerT_bridging::SetEquationSystem",
//			"group %d has external force so equations cannot be reset", group+1);
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

/* initiate the process of writing output from all output sets */
void FEManagerT_bridging::WriteOutput(double time)
{
	/* modified element status */
	if (fElementStatus.Length() > 0)
	{
		/* the continuum element solving the coarse scale */
		const ContinuumElementT* coarse = fDrivenCellData.ContinuumElement();
		ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
		if (non_const_coarse)
		{
			/* enable all coarse scale elements for output */
			ArrayT<ElementCardT::StatusT> status(non_const_coarse->NumElements());
			status = ElementCardT::kON;
			non_const_coarse->SetStatus(status);
		}

		/* inherited - write output */
		FEManagerT::WriteOutput(time);

		/* restore active element map */
		if (non_const_coarse) non_const_coarse->SetStatus(fElementStatus);
	}
	else /* inherited - write output */
		FEManagerT::WriteOutput(time);
}

/* write results for a single output set */
void FEManagerT_bridging::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const
{
	/* inherited */
	FEManagerT::WriteOutput(ID, n_values, e_values);
}

/* write a snapshot */
void FEManagerT_bridging::WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
	const dArray2DT& values, const ArrayT<StringT>& labels) const
{
	/* inherited */
	FEManagerT::WriteOutput(file, coords, node_map, values, labels);
}

/* initialize the ghost node information */
void FEManagerT_bridging::InitGhostNodes(const StringT& field, const ArrayT<StringT>& ghost_id_list, bool include_image_nodes)
{
	const char caller[] = "FEManagerT_bridging::InitGhostNodes";

	/* collect ghost nodes */
	fModelManager->ManyNodeSets(ghost_id_list, fGhostNodes);

	/* get atomistic field */
	FieldT* the_field = fNodeManager->Field(field);
	if (!the_field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", field.Pointer());

	/* create controller to driver solution */
	if (!fSolutionDriver) {
	
		/* construct new contoller */
		fSolutionDriver = fNodeManager->NewKBC_Controller(*the_field, KBC_ControllerT::kPrescribed);

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
			KBC_cards[dex++].SetValues(fGhostNodes[i], j, KBC_CardT::kDsp, NULL, 0.0);

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
	if (fLogging == GlobalT::kVerbose) {
		fGhostNodes++;
		fMainOut << "\n Ghost nodes: " << fGhostNodes.Length() << '\n';
		fMainOut << fGhostNodes.wrap(5) << endl;
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
		is_ghost = 0;

//NOTE: for atom decomp, reproducing crystal everywhere means we cannot
//      leave atoms not owned by this processor marked with 1
#if 0
		/* initialize potential non-ghost nodes */		
		is_ghost = 1;
		const int* p = part_nodes->Pointer();
		int npn = part_nodes->Length();
		for (int i = 0; i < npn; i++)
			is_ghost[*p++] = 0;
#endif
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
		int dof = card.DOF();
		const ScheduleT* schd = card.Schedule();
	
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
void FEManagerT_bridging::InitInterpolation(const StringT& field, const iArrayT& nodes,
	const dArray2DT& coordinates)
{
#pragma unused(field)

	fMainOut << "\n Number of interpolation points. . . . . . . . . = " << nodes.Length() << '\n';

	/* compute interpolation data */
	BridgingScale().InitInterpolation(nodes, &coordinates, NULL, fFollowerCellData);
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
void FEManagerT_bridging::InitProjection(const StringT& field, CommManagerT& comm, const iArrayT& nodes, 
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
		fSolutionDriver = fNodeManager->NewKBC_Controller(*the_field, KBC_ControllerT::kPrescribed);

		/* add to field */
		the_field->AddKBCController(fSolutionDriver);
	}

	/* collect list of projected nodes */
	BridgingScale().CollectProjectedNodes(fDrivenCellData, fProjectedNodes);

	/* generate KBC cards - all degrees of freedom */
	int ndof = the_field->NumDOF();
	ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
	if (make_inactive)
	{
		KBC_cards.Dimension(fProjectedNodes.Length()*ndof);
		int dex = 0;
		for (int j = 0; j < ndof; j++)
			for (int i = 0; i < fProjectedNodes.Length(); i++)
				KBC_cards[dex++].SetValues(fProjectedNodes[i], j, KBC_CardT::kDsp, NULL, 0.0);
	}
	else
		KBC_cards.Dimension(0);

	/* dimension work space */
	fProjection.Dimension(fProjectedNodes.Length(), ndof);
	
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
void FEManagerT_bridging::Project(const dArray2DT& fine_values, dArray2DT& nodal_values)
{
	const char caller[] = "FEManagerT_bridging::Project";

	/* should do some dimension checking */
	BridgingScale().ProjectField(fDrivenCellData, fine_values, nodal_values);
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
	SetFieldValues(field, fProjectedNodes, order, fProjection);
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
	SetFieldValues(field, fProjectedNodes, order, fProjection);
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
  NLSolver* solver = TB_DYNAMIC_CAST(NLSolver*, fSolvers[group]);

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

void FEManagerT_bridging::CollectOverlapRegion_free(iArrayT& overlap_cell, iArrayT& overlap_node) const
{
	const char caller[] = "FEManagerT_bridging::CollectOverlapRegion_free";

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
	for (int i = 0; i < fProjectedNodes.Length(); i++) 
		is_overlap_node[fProjectedNodes[i]] = 'f';

	/* find cells in overlap region */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();
	ArrayT<char> is_overlap_cell(nel);
	is_overlap_cell = 'f';
	int num_overlap_cell = 0;

	/* should really only include cells which contains follower points bonded to
	 * free points. This includes any cells containing followers. */
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	for (int i = 0; i < interpolating_cell.Length(); i++) {
		int cell = interpolating_cell[i];
		if (is_overlap_cell[cell] == 'f')
			num_overlap_cell++;
		is_overlap_cell[cell] = 't';
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
			
	if (fLogging == GlobalT::kVerbose) {
		overlap_cell++;
		fMainOut << "\n overlap cells: " << overlap_cell.Length() << '\n';
		fMainOut << overlap_cell.wrap(5) << endl;
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
	
	for (int i = 0; i < fProjectedNodes.Length(); i++) /* must be a free node */ {
		char& t_f = is_overlap_node[fProjectedNodes[i]];
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

	if (fLogging == GlobalT::kVerbose) {
		overlap_node++;
		fMainOut << "\n overlap nodes: " << overlap_node.Length() << '\n';
		fMainOut << overlap_node.wrap(5) << endl;
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
		ParticlePairT* particle_pair = dynamic_cast<ParticlePairT*>(element_base);
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

/* compute contribution from bonds to ghost atoms */
void FEManagerT_bridging::ComputeSum_signR_Na(const dArrayT& R_i, const RaggedArray2DT<int>& ghost_neighbors, 
	const dArray2DT& coords, const InverseMapT& overlap_node_map, dArrayT& sum_R_N,
	AutoArrayT<int>& overlap_cell_i) const
{
	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail("FEManagerT_bridging::ComputeSum_signR_Na", "interpolation data not set");

	/* interpolating data */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const dArray2DT& interpolating_weights = fFollowerCellData.InterpolationWeights();

#ifdef __DEBUG__
	AutoArrayT<int> a0;
	AutoArrayT<int> a1;
#endif

	/* accumulate ghost bond contribution */
	sum_R_N = 0.0;
	overlap_cell_i.Dimension(0);
	double R = R_i.Magnitude();
	dArrayT bond(R_i.Length());
	iArrayT neighbors;
	dArrayT weights;
	for (int i = 0; i < ghost_neighbors.MajorDim(); i++)
		if (ghost_neighbors.MinorDim(i) > 1) /* not just self */ {
		
			/* get neighbors */
			ghost_neighbors.RowAlias(i, neighbors);
			
			/* search neighbors for bonds matching R_i */
			for (int j = 1; j < neighbors.Length(); j++) {
			
				/* bond terminating at follower node */
				bond.DiffOf(coords(neighbors[j]), coords(neighbors[0]));
				double L_b = bond.Magnitude();
			
				/* bond direction and length */
				double cosRR = dArrayT::Dot(R_i, bond)/L_b/R;
				if (fabs(R - L_b)/R < 1.0e-02 && fabs(fabs(cosRR) - 1.0) < kSmall) {

					/* cell containing the point */
					int follower_point_index = follower_point_map.Map(neighbors[j]);
					int cell = interpolating_cell[follower_point_index];

					/* interpolating weights */
					interpolating_weights.RowAlias(follower_point_index, weights);
					const iArrayT& nodes_U = coarse->ElementCard(cell).NodesU();

					/* accumulate contributions to free nodes */
					bool has_overlap_node = false;
					for (int k = 0; k < nodes_U.Length(); k++) {
						int overlap_node_index = overlap_node_map.Map(nodes_U[k]);
						if (overlap_node_index > -1) {
							sum_R_N[overlap_node_index] += cosRR*weights[k]*L_b/R;
							has_overlap_node = true;
						}
					}
					
					/* collect cells */
					if (has_overlap_node) {
						overlap_cell_i.AppendUnique(cell);

#ifdef __DEBUG__
						/* collect bond atoms */
						a0.Append(neighbors[0]);
						a1.Append(neighbors[j]);
#endif
					}
				}			
			}
		}

#ifdef __DEBUG__
	fMainOut << "\n contributing bonds:\n";
	for (int i = 0; i < a0.Length(); i++)
		fMainOut << a0[i]+1 << " " << a1[i]+1 << '\n';
	fMainOut.flush();
#endif
}

/* compute contribution from bonds to ghost atoms */
void FEManagerT_bridging::ComputeSum_signR_Na(const dArrayT& R_i, const RaggedArray2DT<int>& ghost_neighbors, 
	const dArray2DT& coords, const InverseMapT& overlap_node_map, dArrayT& sum_R_N) const
{
	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail("FEManagerT_bridging::ComputeSum_signR_Na", "interpolation data not set");

	/* interpolating data */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const dArray2DT& interpolating_weights = fFollowerCellData.InterpolationWeights();

	/* accumulate ghost bond contribution */
	sum_R_N = 0.0;
	double R = R_i.Magnitude();
	dArrayT bond(R_i.Length());
	iArrayT neighbors;
	dArrayT weights;
	for (int i = 0; i < ghost_neighbors.MajorDim(); i++)
		if (ghost_neighbors.MinorDim(i) > 1) /* not just self */ {
		
			/* get neighbors */
			ghost_neighbors.RowAlias(i, neighbors);
			
			/* search neighbors for bonds matching R_i */
			for (int j = 1; j < neighbors.Length(); j++) {
			
				/* bond terminating at follower node */
				bond.DiffOf(coords(neighbors[j]), coords(neighbors[0]));
				double L_b = bond.Magnitude();
			
				/* bond direction and length */
				double cosRR = dArrayT::Dot(R_i, bond)/L_b/R;
				if (fabs(R - L_b)/R < 1.0e-02 && fabs(fabs(cosRR) - 1.0) < kSmall) {

					/* cell containing the point */
					int follower_point_index = follower_point_map.Map(neighbors[j]);
					int cell = interpolating_cell[follower_point_index];

					/* interpolating weights */
					interpolating_weights.RowAlias(follower_point_index, weights);
					const iArrayT& nodes_U = coarse->ElementCard(cell).NodesU();

					/* accumulate contributions to free nodes */
					for (int k = 0; k < nodes_U.Length(); k++) {
						int overlap_node_index = overlap_node_map.Map(nodes_U[k]);
						if (overlap_node_index > -1)
							sum_R_N[overlap_node_index] += cosRR*weights[k]*L_b/R;
					}
				}		
			}
		}
}

void FEManagerT_bridging::Compute_df_dp(const dArrayT& R, double V_0,
	const ArrayT<char>& cell_type, const InverseMapT& overlap_cell_map, const InverseMapT& overlap_node_map, 
	const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, LAdMatrixT& ddf_dpdp) const
{
	const char caller[] = "FEManagerT_bridging::Compute_df_dp_1";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");

	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse->NumElementNodes();
	int nip = coarse->NumIP();

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse->ShapeFunction(), element_coords);
	shapes.Initialize();

	/* integrate bond density term */
	dArrayT rho_1(nip);
	rho_1 = 1.0;
	dMatrixT grad_Na(nsd,nen);
	for (int i = 0; i < cell_type.Length(); i++)
		if (cell_type[i] != p_0) /* non-zero bond density */ 
		{
			/* set element information */
			const iArrayT& nodesX = coarse->ElementCard(i).NodesX();
			element_coords.SetLocal(nodesX);
			shapes.SetDerivatives();

			/* bond density */
			const double* p = NULL;
			if (cell_type[i] == p_x) {
				int overlap_cell_index = overlap_cell_map.Map(i);
				if (overlap_cell_index == -1) ExceptionT::GeneralFail(caller);
				p = rho(overlap_cell_index);
			}
			else /* rho = 1.0 */
				p = rho_1.Pointer();

			/* integration parameters */
			const double* j = shapes.IPDets();
			const double* w = shapes.IPWeights();
		
			/* integrate */
			shapes.TopIP();
			while (shapes.NextIP()) {

				/* get shape function gradients */
				shapes.GradNa(grad_Na);
		
				/* integration factor */
				double jw = (*j++)*(*w++);
				double p_jw_by_V = (*p++)*jw/V_0;
		
				/* loop over nodes */
				for (int k = 0; k < nodesX.Length(); k++) {
			
					int overlap_node_index = overlap_node_map.Map(nodesX[k]);
					if (overlap_node_index > -1) /* node is in overlap */ {
				
						/* inner product of bond and shape function gradient */
						double R_dot_dN = grad_Na.DotCol(k, R);
				
						/* assemble */
						f_a[overlap_node_index] += (R_dot_dN*p_jw_by_V);
					}
				}
			}
		}

	/* gradient work space */
	const ParentDomainT& parent_domain = shapes.ParentDomain();	
	ArrayT<dMatrixT> ip_gradient(nip);
	dMatrixT jacobian_inv(nsd);
	dMatrixT A(nsd,nip), ATA(nip), ATA_int(nip);
	for (int i = 0; i < nip; i++) {
		ip_gradient[i].Dimension(nsd,nip);
		parent_domain.IPGradientTransform(i, ip_gradient[i]);
	}

	//TEMP - coupling all integration points
	dArray2DT df_a_dp(f_a.Length(), df_dp.Length()); // BIG array !!!!!!!
	df_a_dp = 0.0;
	dArray2DT df_dp_smoothing = df_dp;
	df_dp_smoothing = 0.0;

	/* compute residual */
	dArray2DT df(nen, nip);
	dMatrixT ddp_i_dpdp(nip);
	df_dp = 0.0;
	ddf_dpdp = 0.0;
	dArrayT element_rho;
	dArrayT element_force;
	for (int i = 0; i < cell_type.Length(); i++)
		if (cell_type[i] == p_x) /* unknown bond density */ 
		{
			/* index within list of overlap cells */
			int overlap_cell_index = overlap_cell_map.Map(i);
			if (overlap_cell_index == -1) ExceptionT::GeneralFail(caller);
		
			/* set element information */
			const iArrayT& nodesX = coarse->ElementCard(i).NodesX();
			element_coords.SetLocal(nodesX);
			shapes.SetDerivatives();

			/* integration parameters */
			const double* p = rho(overlap_cell_index);
			const double* j = shapes.IPDets();
			const double* w = shapes.IPWeights();
		
			/* integrate */
			ATA_int = 0.0;
			ddp_i_dpdp = 0.0;
			df = 0.0;
			shapes.TopIP();
			while (shapes.NextIP()) {

				/* integration factor */
				int ip = shapes.CurrIP();
				double jw = (*j++)*(*w++);
				double pm1_jw = ((*p++) - 1)*jw; /* penalization */
				double jw_by_V = jw/V_0;
	
				/* get shape function gradients */
				shapes.GradNa(grad_Na);

				/* integrate density gradient matrix */
				parent_domain.DomainJacobian(element_coords, ip, jacobian_inv);
				jacobian_inv.Inverse();
				A.MultATB(jacobian_inv, ip_gradient[ip]);
				ATA.MultATB(A,A);
				ATA_int.AddScaled(smoothing*jw, ATA);
			
				/* add penalized force term */
				df_dp(overlap_cell_index,ip) += k2*pm1_jw;
				ddp_i_dpdp(ip,ip) += k2*jw;
					
				/* integrate the bond density term over the element */
				for (int k = 0; k < nodesX.Length(); k++) {
			
					int overlap_node_index = overlap_node_map.Map(nodesX[k]);
					if (overlap_node_index > -1) /* node is in overlap */ {
				
						/* inner product of bond and shape function gradient */
						double R_dot_dN = grad_Na.DotCol(k, R);
				
						/* assemble */
						df(k, ip) += (R_dot_dN*jw_by_V);
					}
				}
			}
		
			/* accumulate */
			for (int j = 0; j < nodesX.Length(); j++) {
			
				int overlap_node_index = overlap_node_map.Map(nodesX[j]);
				if (overlap_node_index > -1) /* node is in overlap */ {
			
					/* across all integration points */
					for (int k = 0; k < nip; k++) {
				
						/* force */
						df_dp(overlap_cell_index,k) += f_a[overlap_node_index]*df(j,k);

						//TEMP - accumulate
						df_a_dp(overlap_node_index, overlap_cell_index*nip + k) += df(j,k);
					}
				}
			}
		
			/* regularization contribution to force and stiffness */
			df_dp.RowAlias(overlap_cell_index, element_force);
			rho.RowAlias(overlap_cell_index, element_rho);		
			ATA_int.Multx(element_rho, element_force, 1.0, dMatrixT::kAccumulate);

			/* penalty regularization */
			ATA_int += ddp_i_dpdp;
		
			ddf_dpdp.AddBlock(overlap_cell_index*nip, overlap_cell_index*nip, ATA_int);
		}

	/* add coupling term */
	for (int i = 0; i < df_a_dp.MajorDim(); i++)
		ddf_dpdp.Outer(df_a_dp(i), df_a_dp(i), 1.0, dMatrixT::kAccumulate);
}

/* compute reduced connectivity list */
void FEManagerT_bridging::GhostNodeBonds(const RaggedArray2DT<int>& neighbors, RaggedArray2DT<int>& ghost_neighbors, 
	InverseMapT& overlap_cell_map) const
{
	/* cell containing each interpolating point */
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();

	/* keep only bonds to follower points */
	iArrayT neighbors_i;
	AutoFill2DT<int> gh_neigh(neighbors.MajorDim(), 1, 25, 10);
	for (int i = 0; i < neighbors.MajorDim(); i++) {
		neighbors.RowAlias(i, neighbors_i);
		
		/* self is first neighbor */
		gh_neigh.Append(i, neighbors_i[0]);
			
		/* search other neighbors */
		for (int j = 1; j < neighbors_i.Length(); j++) {
		
			/* the neighbor */
			int neighbor = neighbors_i[j];
			
			/* point is follower */
			int neighbor_local = follower_point_map.Map(neighbor);
			if (neighbor_local > -1) {
			
				/* cell containing the point */
				int cell = interpolating_cell[neighbor_local];

				/* cell is in overlap region */
				if (overlap_cell_map.Map(cell) > -1)
					gh_neigh.Append(i, neighbor);
			}
		}
	}

	/* copy/compress */
	ghost_neighbors.Copy(gh_neigh);
}

void FEManagerT_bridging::GhostNodeBonds(const RaggedArray2DT<int>& neighbors,
	RaggedArray2DT<int>& ghost_neighbors, iArrayT& overlap_cell) const
{
	const char caller[] = "FEManagerT_bridging::GhostNodeBonds";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");
	int nel = coarse->NumElements();
	ArrayT<char> is_overlap_cell(nel);
	is_overlap_cell = 'f';

	/* interpolation data */
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();

	/* keep only bonds (from free atoms) to ghost points */
	int num_overlap = 0;
	iArrayT neighbors_i;
	AutoFill2DT<int> gh_neigh(neighbors.MajorDim(), 1, 25, 10);
	for (int i = 0; i < neighbors.MajorDim(); i++) {

		/* the neighbor list */
		neighbors.RowAlias(i, neighbors_i);
		
		/* self is first neighbor */
		gh_neigh.Append(i, neighbors_i[0]);
			
		/* search other neighbors for ghosts */
		for (int j = 1; j < neighbors_i.Length(); j++) {
		
			/* the neighbor */
			int neighbor = neighbors_i[j];
			
			/* point is follower */
			int neighbor_local = follower_point_map.Map(neighbor);
			if (neighbor_local > -1) {
			
				/* keep neighbor */
				gh_neigh.Append(i, neighbor);
			
				/* mark cell */
				int cell = interpolating_cell[neighbor_local];
				if (cell < 0) ExceptionT::GeneralFail(caller, "no cell for point %d", neighbor+1);
				if (is_overlap_cell[cell] == 'f') num_overlap++;
				is_overlap_cell[cell] = 't';
			}
		}
	}

	/* copy/compress */
	ghost_neighbors.Copy(gh_neigh);
	
	/* collect overlap cells */
	overlap_cell.Dimension(num_overlap);
	num_overlap = 0;
	for (int i = 0; i < is_overlap_cell.Length(); i++)
		if (is_overlap_cell[i] == 't')
			overlap_cell[num_overlap++] = i;

//NOTE: up to this point, overlap cells are classified as those containing an active ghost
//      atom bond; while overlap nodes are those nodes whose support includes a ghost atom.
//      Now, we grow the list of overlap cells to include the entire support of all nodes
//      whose support contains an active ghost node bond. This must be consistent with the
//      other version of FEManagerT_bridging::GhostNodeBonds used to collect the entire,
//      potential overlap region.
#if 0
	/* mark overlap nodes - nodes whose support contains an active
	 * ghost atom bond */
	int nnd = fNodeManager->NumNodes();	
	ArrayT<char> is_overlap_node(nnd);
	is_overlap_node = 'f';
	num_overlap = 0;
	for (int i = 0; i < overlap_cell.Length(); i++) {
		const iArrayT& nodesX = coarse->ElementCard(overlap_cell[i]).NodesX();
		for (int j = 0; j < nodesX.Length(); j++) {
			int nd = nodesX[j];
			if (is_overlap_node[nd] == 'f') {
				num_overlap++;
				is_overlap_node[nd] = 't';
			}
		}
	}

	/* loop over all elements looking for overlap nodes */
	is_overlap_cell = 'f';
	num_overlap = 0;
	for (int i = 0; i < nel; i++) {
		const iArrayT& nodesX = coarse->ElementCard(i).NodesX();
		bool has_overlap_node = false;
		for (int j = 0; !has_overlap_node && j < nodesX.Length(); j++)
			if (is_overlap_node[nodesX[j]] == 't')
				has_overlap_node = true;
	
		if (has_overlap_node) {
			is_overlap_cell[i] = 't';
			num_overlap++;
		}
	}

	/* collect overlap cells */
	overlap_cell.Dimension(num_overlap);
	num_overlap = 0;
	for (int i = 0; i < is_overlap_cell.Length(); i++)
		if (is_overlap_cell[i] == 't')
			overlap_cell[num_overlap++] = i;
#endif
}

/* collect list of cells not containing any active bonds */
void FEManagerT_bridging::BondFreeElements(const RaggedArray2DT<int>& ghost_neighbors_i, 
	AutoArrayT<int>& bondfree_cell_i) const
{
	const char caller[] = "FEManagerT_bridging::BondFreeElements";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");
	int nel = coarse->NumElements();
	ArrayT<char> is_bond_free(nel);
	is_bond_free = 't';

	/* mark cells containing projection points */
	int num_bond_free = nel;
	const RaggedArray2DT<int>& point_in_cell = fDrivenCellData.PointInCell();
	for (int i = 0; i < nel; i++)
		if (point_in_cell.MinorDim(i) > 0) {
			is_bond_free[i] = 'f';
			num_bond_free--;
		}

	/* interpolation data */
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();

	/* process bonds to ghost points */
	iArrayT neighbors_i;
	for (int i = 0; i < ghost_neighbors_i.MajorDim(); i++) {

		/* the neighbor list */
		ghost_neighbors_i.RowAlias(i, neighbors_i);
			
		/* search other neighbors for ghosts */
		for (int j = 1; j < neighbors_i.Length(); j++) {

				/* mark cell */
				int neighbor_local = follower_point_map.Map(neighbors_i[j]);	
				int cell = interpolating_cell[neighbor_local];
				if (is_bond_free[cell] == 't') {
					is_bond_free[cell] = 'f';
					num_bond_free--;
				}
			}
		}

	/* collect bond-free cells */
	bondfree_cell_i.Dimension(num_bond_free);
	num_bond_free = 0;
	for (int i = 0; i < is_bond_free.Length(); i++)
		if (is_bond_free[i] == 't')
			bondfree_cell_i[num_bond_free++] = i;
}

/* count number of bonds terminating with the domain of each element integration point */
void FEManagerT_bridging::CountIPBonds(const ArrayT<int>& elements, iArray2DT& ip_counts) const
{
	const char caller[] = "FEManagerT_bridging::CountIPBonds";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");

	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse->NumElementNodes();
	int nip = ip_counts.MinorDim();

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse->GeometryCode(), nip, element_coords);
	shapes.Initialize();
	const ParentDomainT& parent_domain = shapes.ParentDomain();

	/* interpolation information */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const RaggedArray2DT<double>& point_in_cell_coords = fFollowerCellData.PointInCellCoords();
	int nel = point_in_cell.MajorDim();

	/* initialize */
	ip_counts = 0;
	dArray2DT point_coords;
	dArrayT point;
	for (int i = 0; i < elements.Length(); i++) {
		int e = elements[i];
		if (point_in_cell.MinorDim(e) > 0)
		{
			point_coords.Alias(point_in_cell.MinorDim(e), nsd, point_in_cell_coords(e));
			for (int j = 0; j < point_coords.MajorDim(); j++)
			{
				point_coords.RowAlias(j, point);
				int ip = parent_domain.IPDomain(point);
				ip_counts(i,ip)++;
			}
		}
	}

	/* verbose */
	if (fLogging == GlobalT::kVerbose) {
		fMainOut << "\n Count of bonds in integration point domains:\n";
		fMainOut << ip_counts << '\n';
		fMainOut.flush();
	}
}

/* generate "inverse" connectivities for active elements in the support of active nodes */
void FEManagerT_bridging::TransposeConnects(const ContinuumElementT& element_group, 
	const ArrayT<int>& active_nodes, const ArrayT<int>& active_elements, 
	RaggedArray2DT<int>& transpose_connects)
{
	InverseMapT active_node_map;
	active_node_map.SetOutOfRange(InverseMapT::MinusOne);
	active_node_map.SetMap(active_nodes);

	/* generate "inverse" connectivities */	
	AutoFill2DT<int> elements_per_node(active_nodes.Length(), 1, 25, 10);
	for (int i= 0; i < active_elements.Length(); i++) {
		const iArrayT& nodes = element_group.ElementCard(active_elements[i]).NodesU();
		for (int j = 0; j < nodes.Length(); j++) {
			int active_node_index = active_node_map.Map(nodes[j]);
			if (active_node_index != -1)
				elements_per_node.Append(active_node_index, active_elements[i]);
		}
	}
	
	/* copy into return value */
	transpose_connects.Copy(elements_per_node);
}

/* group bonds into shells */
int FEManagerT_bridging::NumberShells(const dArray2DT& bonds, iArrayT& shell, dArrayT& shell_bond_length) const
{
	/* initialize shell numbering */
	shell.Dimension(bonds.MajorDim());
	shell = -1;

	/* run through bonds and group */
	AutoArrayT<double> bond_length;
	dArrayT bond;
	int shell_count = 0;
	for (int i = 0; i < bonds.MajorDim(); i++)
		if (shell[i] == -1)
		{
			/* new shell */
			shell[i] = shell_count;
			bonds.RowAlias(i, bond);
			double r = bond.Magnitude();
			bond_length.Append(r);
			
			/* look for bonds with same length */
			for (int j = i+1; j < bonds.MajorDim(); j++)
			{
				bonds.RowAlias(j, bond);
				if (fabs(r - bond.Magnitude()) < kSmall)
					shell[j] = shell_count;
			}

			shell_count++;
		}

	bond_length.Swap(shell_bond_length);
	return shell_count;
}

#endif  /* BRIDGING_ELEMENT */
