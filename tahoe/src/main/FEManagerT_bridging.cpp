/* $Id: FEManagerT_bridging.cpp,v 1.16.4.17 2004-04-24 01:43:37 paklein Exp $ */
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
#include "Chain1D.h"
#include "BondLatticeT.h"
#include "ShapeFunctionT.h"
#include "nArrayGroupT.h"
#include "nVariMatrixT.h"
//TEMP
#include "LAdMatrixT.h"
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

void FEManagerT_bridging::CorrectOverlap(const RaggedArray2DT<int>& point_neighbors, const dArray2DT& point_coords, double smoothing, double k2)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap";

	/* finding free vs projected nodes */
	int nnd = fNodeManager->NumNodes();
	ArrayT<char> node_type(nnd);
	node_type = free_;
	for (int i = 0; i < fProjectedNodes.Length(); i++) /* mark projected nodes */
		node_type[fProjectedNodes[i]] = not_free_;

	/* map describing free or ghost points */
	ArrayT<char> point_type(point_coords.MajorDim());
	point_type = free_;
	const ArrayT<int>& point_in_cell_data = fFollowerCellData.PointInCell().Data();
	for (int i = 0; i < point_in_cell_data.Length(); i++)
		point_type[point_in_cell_data[i]] = not_free_;

	/* collect only bonds terminating with ghost points */
	RaggedArray2DT<int> ghost_neighbors_all;
	iArrayT overlap_cell_all;
	InverseMapT overlap_cell_all_map;
	overlap_cell_all_map.SetOutOfRange(InverseMapT::MinusOne);	
	GhostNodeBonds(point_neighbors, ghost_neighbors_all, overlap_cell_all);
	overlap_cell_all_map.SetMap(overlap_cell_all);
	if (fPrintInput) {
		fMainOut << "\n Bonds to interpolation points (self as leading neighbor):\n";
		fMainOut << setw(kIntWidth) << "row" << "  n..." << '\n';
		iArrayT tmp(ghost_neighbors_all.Length(), ghost_neighbors_all.Pointer());
		tmp++;
		ghost_neighbors_all.WriteNumbered(fMainOut);
		tmp--;
		fMainOut.flush();
	}

	/* coarse scale element group */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	int nel = coarse->NumElements();

	/* Cauchy-Born constitutive model */
	const MaterialListT& mat_list = coarse->MaterialsList();
	if (mat_list.Length() > 1) ExceptionT::GeneralFail(caller, "expecting only 1 material not %d", mat_list.Length());
	ContinuumMaterialT* cont_mat = mat_list[0];
	Hex2D* hex_2D = dynamic_cast<Hex2D*>(cont_mat);
	FCC3D* fcc_3D = dynamic_cast<FCC3D*>(cont_mat);
	Chain1D* mat_1D = dynamic_cast<Chain1D*>(cont_mat);
	if (!mat_1D && !hex_2D && !fcc_3D) ExceptionT::GeneralFail(caller, "could not resolve C-B material");

	/* lattice information */
	const BondLatticeT& bond_lattice = (hex_2D) ? hex_2D->BondLattice() : ((fcc_3D) ? fcc_3D->BondLattice() : mat_1D->BondLattice());
	const dArray2DT& bonds = bond_lattice.Bonds();
	double V_0 = (hex_2D) ? hex_2D->CellVolume() : ((fcc_3D) ? fcc_3D->CellVolume() : mat_1D->CellVolume());
	double R_0 = (hex_2D) ? hex_2D->NearestNeighbor() : ((fcc_3D) ? fcc_3D->NearestNeighbor() : mat_1D->NearestNeighbor());

	/* unknown bond densities */
	dArray2DT bond_densities(overlap_cell_all.Length(), coarse->NumIP()*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space that changes for each bond family */
	LAdMatrixT ddf_dpdp_i;
	CCSMatrixT ddf_dpdp_i_(Output(), GlobalMatrixT::kZeroPivots);
	
	nVariMatrixT<double> ddf_dpdp_i_man(0, ddf_dpdp_i);
	dArray2DT p_i, dp_i, df_dp_i;
	nArray2DGroupT<double> ip_unknown_group(0, false, coarse->NumIP());
	ip_unknown_group.Register(p_i);
	ip_unknown_group.Register(dp_i);
	ip_unknown_group.Register(df_dp_i);
	dArrayT sum_R_N;
	dArrayT f_a;
	nArrayGroupT<double> overlap_node_group(0, false);
	overlap_node_group.Register(sum_R_N);
	overlap_node_group.Register(f_a);

	/* solve unknowns */
	ArrayT<char> cell_type(nel);
	dArrayT R_i(point_coords.MinorDim());
	AutoArrayT<int> overlap_cell_i;	
	InverseMapT overlap_cell_i_map;
	overlap_cell_i_map.SetOutOfRange(InverseMapT::MinusOne);
	AutoArrayT<int> overlap_node_i;
	InverseMapT overlap_node_i_map;
	overlap_node_i_map.SetOutOfRange(InverseMapT::MinusOne);
	iArray2DT bond_densities_i_eq;
	nVariArray2DT<int> bond_densities_i_eq_man(0, bond_densities_i_eq, coarse->NumIP());
	RaggedArray2DT<int> inv_connects_i; 
	RaggedArray2DT<int> inv_equations_i; 
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond family at a time */ {

		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
	
		/* collect "acive" ghost node bonds - bonds of type R_i terminating at a ghost node */
		RaggedArray2DT<int> ghost_neighbors_i;
		GhostNodeBonds(R_i, point_coords, ghost_neighbors_all, ghost_neighbors_i, overlap_cell_i, overlap_node_i);
		overlap_cell_i_map.SetMap(overlap_cell_i);
		overlap_node_i_map.SetMap(overlap_node_i);
		if (fPrintInput) {
			fMainOut << "overlap cells for bond: " << ": {" << R_i.no_wrap() << "}:\n";
			iArrayT tmp;
			tmp.Alias(overlap_cell_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;
		}

		/* classify element types */
		cell_type = p_1; /* all full density */
		for (int j = 0; j < nel; j++) /* no density */ {
			const iArrayT& nodes = coarse->ElementCard(j).NodesX();
			bool no_density = true;
			for (int k = 0; no_density && k < nodes.Length(); k++)
				if (node_type[nodes[k]] == free_)
					no_density = false;
			
			/* no free nodes */
			if (no_density)
				cell_type[j] = p_0;
		}
		for (int j = 0; j < overlap_cell_i.Length(); j++) /* unknown density */
			cell_type[overlap_cell_i[j]] = p_x;
			
		/* dimension work space */
		ddf_dpdp_i_man.SetDimensions(overlap_cell_i.Length()*coarse->NumIP());
		ip_unknown_group.SetMajorDimension(overlap_cell_i.Length(), false);
		overlap_node_group.Dimension(overlap_node_i.Length(), false);
		
		/* number unknown bond densities */
		bond_densities_i_eq_man.SetMajorDimension(overlap_cell_i.Length(), false);
		bond_densities_i_eq = -1;
		int num_eq = 0;
		for (int j = 0; j < overlap_cell_i.Length(); j++)
			for (int k = 0; k < bond_densities_i_eq.MinorDim(); k++)
				bond_densities_i_eq(j,k) = ++num_eq;

		/* "inverse" connectivities for overlap cells in the support of overlap nodes */
		TransposeConnects(*coarse, overlap_node_i, overlap_cell_i, inv_connects_i);
		
		/* collect equation numbers for "inverse" elements */
		inv_equations_i.Configure(inv_connects_i, coarse->NumIP());
		for (int j = 0; j < inv_connects_i.MajorDim(); j++) {
			const int* element_per_node = inv_connects_i(j);
			int* equations_per_node = inv_equations_i(j);
			for (int k = 0; k < inv_connects_i.MinorDim(j); k++) {
				int overlap_cell_index = overlap_cell_i_map.Map(element_per_node[k]);
				int* equations = bond_densities_i_eq(overlap_cell_index);
				for (int l = 0; l < bond_densities_i_eq.MinorDim(); l++)
					*equations_per_node++ = *equations++;
			}
		}

		/* configure linear solver */
		ddf_dpdp_i_.AddEquationSet(inv_equations_i);
		ddf_dpdp_i_.Initialize(num_eq, num_eq, 1);
		
		/* compute contribution from bonds terminating at "ghost" atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors_i, point_coords, overlap_node_i_map, sum_R_N);

		/* initialize */
		p_i = 1.0;
		
		/* compute residual - add Cauchy-Born contribution */
		f_a = sum_R_N;
		dArray2DT df_dp_i_(df_dp_i);
		Compute_df_dp(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
			bond_densities_i_eq, inv_connects_i, inv_equations_i,
			p_i, f_a, smoothing, k2, df_dp_i_, ddf_dpdp_i_);

		f_a = sum_R_N;
		Compute_df_dp(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);

#if __DEBUG__
int prec = fMainOut.precision();
fMainOut.precision(12);
fMainOut << "f_old = \n" << df_dp_i << '\n';
fMainOut << "K_old = \n" << ddf_dpdp_i << '\n';
fMainOut << "f_new = \n" << df_dp_i_ << '\n';
fMainOut << "K_new = \n" << ddf_dpdp_i_ << '\n';
fMainOut.precision(prec);
fMainOut.flush();
#endif
		
		/* solve bond densities */
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+06;		
		int max_iter = 5;
		int iter = 0;
		double error_0 = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));		
		double error = error_0;
		while (iter++ < max_iter && error > abs_tol && error/error_0 > rel_tol && error/error_0 < div_tol) {

			/* catch errors in linear solver */
			try {

//TEMP new solver
				dp_i.SetToScaled(-1.0, df_dp_i_);
				dArrayT tmp_;
				tmp_.Alias(dp_i);
				ddf_dpdp_i_.Solve(tmp_);
				fMainOut << "dp_i =\n" << tmp_ << endl;
//TEMP

				/* solve system */
				dp_i.SetToScaled(-1.0, df_dp_i);
				dArrayT tmp;
				tmp.Alias(dp_i);
#if __DEBUG__
				fMainOut << "f:\n" << tmp << endl;
#endif
				ddf_dpdp_i.LinearSolve(tmp);
#if __DEBUG__
				fMainOut << "dp_i =\n" << tmp << endl;
#endif

				/* update densities */
				p_i += dp_i;

				/* recompute residual */			
				f_a = sum_R_N;
				Compute_df_dp(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
				error = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
				cout << setw(kIntWidth) << iter << ": e/e_0 = " << error/error_0 << endl;

			} /* end try */

			catch (ExceptionT::CodeT error) {
				iter = max_iter; /* exit */
			}
		}

		if (error > abs_tol && error/error_0 > rel_tol) /* converged */ {
			cout << "FAIL" << endl;
			p_i = 1.0;
		}

		/* save result */
		int nb = bonds.MajorDim();
		for (int k = 0; k < overlap_cell_i.Length(); k++) {
			int cell = overlap_cell_all_map.Map(overlap_cell_i[k]);
			for (int j = 0; j < coarse->NumIP(); j++)
				bond_densities(cell, i+j*nb) = p_i(k,j);
		}		
	}
	
	/* write densities */
	if (fPrintInput) {
	
		/* dimensions */
		const NodeManagerT* node_manager = NodeManager();
		int nsd = node_manager->NumSD();
		int nip = coarse->NumIP();
		int nen = coarse->NumElementNodes();

		/* shape functions */
		LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
		NodeManager()->RegisterCoordinates(element_coords);
		ShapeFunctionT shapes = ShapeFunctionT(coarse->ShapeFunction(), element_coords);
		shapes.Initialize();
	
		/* collect integration point coordinates */
		dArray2DT coords(overlap_cell_all.Length()*nip, nsd);
		dArrayT ip_coords;
		int index = 0;
		for (int i = 0; i < overlap_cell_all.Length(); i++) {
		
			/* collect element coordinates */
			const iArrayT& nodes = coarse->ElementCard(overlap_cell_all[i]).NodesX();
			element_coords.SetLocal(nodes);
			
			/* compute ip coordinates */
			for (int j = 0; j < nip; j++) {
				coords.RowAlias(index++, ip_coords);
				shapes.IPCoords(ip_coords, j);
			}
		}
		
		/* collect output labels */
		ArrayT<StringT> bond_labels(bonds.MajorDim());
		for (int i = 0; i < bond_labels.Length(); i++)
			bond_labels[i].Append("p_", i+1);
	
		/* output file */
		StringT file;
		file.Root(fMainIn.filename());
		file.Append(".bond_density.out");
			
		/* write output */
		dArray2DT n_values(overlap_cell_all.Length()*nip, bonds.MajorDim(), bond_densities.Pointer()); /* maintain order, change shape */
		iArrayT point_map(coords.MajorDim());
		point_map.SetValueToPosition();
		WriteOutput(file, coords, point_map, n_values, bond_labels);			
	}

	/* write unknowns into the state variable space */
	ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
	for (int i = 0; i < overlap_cell_all.Length(); i++) {
	
		/* element information */
		ElementCardT& element = non_const_coarse->ElementCard(overlap_cell_all[i]);
	
		/* allocate space */
		element.Dimension(0, bond_densities.MinorDim());
		
		/* copy in densities */
		element.DoubleData() = bond_densities(i);
	}	
}

/* compute internal correction for the overlap region */
void FEManagerT_bridging::CorrectOverlap_ghost(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, double smoothing, double k2)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_ghost";

	/* collect nodes and cells in the overlap region */
	iArrayT overlap_cell;
	iArrayT overlap_node;
	CollectOverlapRegion_free(overlap_cell, overlap_node);
	if (overlap_node.Length() == 0) return;

	/* map of overlap_node in local numbering */
	InverseMapT overlap_node_map;
	overlap_node_map.SetOutOfRange(InverseMapT::MinusOne);
	overlap_node_map.SetMap(overlap_node);

	/* map of overlap_cell in local numbering */
	InverseMapT overlap_cell_map;
	overlap_cell_map.SetOutOfRange(InverseMapT::MinusOne);
	overlap_cell_map.SetMap(overlap_cell);

	/* compute reduced connectivity list */
	RaggedArray2DT<int> ghost_neighbors;
	GhostNodeBonds(neighbors, ghost_neighbors, overlap_cell_map);
	if (fPrintInput) {
		fMainOut << "\n Bonds to interpolation points (self as leading neighbor):\n";
		fMainOut << setw(kIntWidth) << "row" << "  n..." << '\n';
		iArrayT tmp(ghost_neighbors.Length(), ghost_neighbors.Pointer());
		tmp++;
		ghost_neighbors.WriteNumbered(fMainOut);
		tmp--;
		fMainOut.flush();
	}

	/* Cauchy-Born constitutive model */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
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
	double V_0 = (hex_2D) ? hex_2D->CellVolume() : ((fcc_3D) ? fcc_3D->CellVolume() : mat_1D->CellVolume());
	double R_0 = (hex_2D) ? hex_2D->NearestNeighbor() : ((fcc_3D) ? fcc_3D->NearestNeighbor() : mat_1D->NearestNeighbor());

	/* unknown bond densities */
	dArray2DT bond_densities(overlap_cell.Length(), coarse->NumIP()*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space */
	//TEMP
	//dArray2DT ddf_dpdp_i(overlap_cell.Length(), coarse->NumIP());
	LAdMatrixT ddf_dpdp_i;
	nVariMatrixT<double> ddf_dpdp_i_man(0, ddf_dpdp_i);
	dArray2DT p_i, dp_i, df_dp_i;
	nArray2DGroupT<double> ip_unknown_group(0, false, coarse->NumIP());
	ip_unknown_group.Register(p_i);
	ip_unknown_group.Register(dp_i);
	ip_unknown_group.Register(df_dp_i);

	/* solve unkowns */
	dArrayT sum_R_N(overlap_node.Length());
	dArrayT f_a(overlap_node.Length());
	dArrayT R_i(coords.MinorDim());
	AutoArrayT<int> overlap_cell_i;
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond density at a time */ {
		
		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
		
		/* compute contribution from bonds to ghost atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors, coords, overlap_node_map, sum_R_N, overlap_cell_i);
		if (fPrintInput) /* debugging */ {
			
			/* coordinates */
			const NodeManagerT* node_manager = NodeManager();
			dArray2DT coords(overlap_node.Length(), node_manager->NumSD());
			coords.RowCollect(overlap_node, node_manager->InitialCoordinates());

			/* output file */
			StringT file;
			file.Root(fMainIn.filename());
			file.Append(".RdN.b", i+1);
			file.Append(".out");
			
			/* write output */
			dArray2DT n_values(overlap_node.Length(), 1, sum_R_N.Pointer());
			ArrayT<StringT> n_labels(1);
			n_labels[0] = "f";
			WriteOutput(file, coords, overlap_node, n_values, n_labels);		
		}
		
		/* dimension the local problem */
		ddf_dpdp_i_man.SetDimensions(overlap_cell_i.Length()*coarse->NumIP());
		ip_unknown_group.SetMajorDimension(overlap_cell_i.Length(), false);

		/* initialize */
		p_i = 1.0;
		
		/* solve for bond densities */
		f_a = sum_R_N;
		Compute_df_dp(R_i, V_0, *coarse, overlap_cell_i, overlap_node_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
		double error_0 = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
		cout << setw(kIntWidth) << i+1 << ": {" << R_i.no_wrap() << "}: " << error_0 << '\n';
		if (fPrintInput) {
			fMainOut << "overlap cells for bond: " << ": {" << R_i.no_wrap() << "}:\n";
			iArrayT tmp;
			tmp.Alias(overlap_cell_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;
		}
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+06;		
		int max_iter = 1000;
		int iter = 0;
		double error = error_0;
		while (iter++ < max_iter && error > abs_tol && error/error_0 > rel_tol && error/error_0 < div_tol) {

			try {		
#if 0
			/* compute update vector (steepest descent) */
			for (int i = 0; i < p_i.Length(); i++) {
				if (fabs(ddf_dpdp_i(i,i)) < kSmall) 
					ExceptionT::BadJacobianDet(caller, "smoothing matrix is singular");
				dp_i[i] = -df_dp_i[i]/ddf_dpdp_i(i,i);
			}
#endif

			dp_i.SetToScaled(-1.0, df_dp_i);
			dArrayT tmp;
			tmp.Alias(dp_i);
#if __DEBUG__
			fMainOut << "f:\n" << tmp << endl;
#endif
			ddf_dpdp_i.LinearSolve(tmp);

			/* update densities */
			p_i += dp_i;
			
			/* recompute residual */			
			f_a = sum_R_N;
			Compute_df_dp(R_i, V_0, *coarse, overlap_cell_i, overlap_node_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
			error = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
			cout << setw(kIntWidth) << iter << ": e/e_0 = " << error/error_0 << endl;
			
			} /* end try */
			catch (ExceptionT::CodeT error) {
				iter = max_iter; /* exit */
			}
		}
		
		if (error > abs_tol && error/error_0 > rel_tol) /* converged */ {
			cout << "FAIL" << endl;
			p_i = 1.0;
		}

		/* save result */
		int nb = bonds.MajorDim();
		for (int k = 0; k < overlap_cell_i.Length(); k++) {
			int cell = overlap_cell_map.Map(overlap_cell_i[k]);
			for (int j = 0; j < coarse->NumIP(); j++)
				bond_densities(cell, i+j*nb) = p_i(k,j);
		}
	}

	/* write densities */
	if (fPrintInput) {
	
		/* dimensions */
		const NodeManagerT* node_manager = NodeManager();
		int nsd = node_manager->NumSD();
		int nip = coarse->NumIP();
		int nen = coarse->NumElementNodes();

		/* shape functions */
		LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
		NodeManager()->RegisterCoordinates(element_coords);
		ShapeFunctionT shapes = ShapeFunctionT(coarse->ShapeFunction(), element_coords);
		shapes.Initialize();
	
		/* collect integration point coordinates */
		dArray2DT coords(overlap_cell.Length()*nip, nsd);
		dArrayT ip_coords;
		int index = 0;
		for (int i = 0; i < overlap_cell.Length(); i++) {
		
			/* collect element coordinates */
			const iArrayT& nodes = coarse->ElementCard(overlap_cell[i]).NodesX();
			element_coords.SetLocal(nodes);
			
			/* compute ip coordinates */
			for (int j = 0; j < nip; j++) {
				coords.RowAlias(index++, ip_coords);
				shapes.IPCoords(ip_coords, j);
			}
		}
		
		/* collect output labels */
		ArrayT<StringT> bond_labels(bonds.MajorDim());
		for (int i = 0; i < bond_labels.Length(); i++)
			bond_labels[i].Append("p_", i+1);
	
		/* output file */
		StringT file;
		file.Root(fMainIn.filename());
		file.Append(".bond_density.out");
			
		/* write output */
		dArray2DT n_values(overlap_cell.Length()*nip, bonds.MajorDim(), bond_densities.Pointer()); /* maintain order, change shape */
		iArrayT point_map(coords.MajorDim());
		point_map.SetValueToPosition();
		WriteOutput(file, coords, point_map, n_values, bond_labels);			
	}

	/* write unknowns into the state variable space */
	ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* element information */
		ElementCardT& element = non_const_coarse->ElementCard(overlap_cell[i]);
	
		/* allocate space */
		element.Dimension(0, bond_densities.MinorDim());
		
		/* copy in densities */
		element.DoubleData() = bond_densities(i);
	}
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

	/* collect list of projected nodes */
	const InterpolationDataT& point_to_node = fDrivenCellData.PointToNode();
	if (point_to_node.Neighbors().MajorDim() > 0) /* check for meshless bridging */ {
		const InverseMapT& driven_node_map = point_to_node.Map();
		driven_node_map.Forward(fProjectedNodes);
	}
	else
		fProjectedNodes.Alias(fDrivenCellData.CellNodes());
	
	/* generate KBC cards - all degrees of freedom */
	int ndof = the_field->NumDOF();
	if (make_inactive)
	{
		ArrayT<KBC_CardT>& KBC_cards = fSolutionDriver->KBC_Cards();
		KBC_cards.Dimension(fProjectedNodes.Length()*ndof);
		int dex = 0;
		for (int j = 0; j < ndof; j++)
			for (int i = 0; i < fProjectedNodes.Length(); i++)
				KBC_cards[dex++].SetValues(fProjectedNodes[i], j, KBC_CardT::kDsp, 0, 0.0);
	}

	/* dimension work space */
	fProjection.Dimension(fProjectedNodes.Length(), ndof);
	
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
			
	if (fPrintInput) {
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

	if (fPrintInput) {
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
			fe->fBridgingScale = dynamic_cast<BridgingScaleT*>(element_base);
		}
		
		/* not found */
		if (!fBridgingScale)
			ExceptionT::GeneralFail("FEManagerT_bridging::BridgingScale",
				"did not find BridgingScaleT element group");
	}
	
	return *fBridgingScale;
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

/* compute Cauchy-Born contribution to the nodal internal force */
void FEManagerT_bridging::Compute_df_dp(const dArrayT& R, double V_0, const ContinuumElementT& coarse, 
	const ArrayT<int>& overlap_cell, const InverseMapT& overlap_node_map, const dArray2DT& rho, 
	dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, LAdMatrixT& ddf_dpdp) const
{
	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse.NumElementNodes();
	int nip = coarse.NumIP();

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse.ShapeFunction(), element_coords);
	shapes.Initialize();

	/* integrate bond density term over cells in overlap */
	dMatrixT grad_Na(nsd,nen);
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* set element information */
		const iArrayT& nodesX = coarse.ElementCard(overlap_cell[i]).NodesX();
		element_coords.SetLocal(nodesX);
		shapes.SetDerivatives();

		/* integration parameters */
		const double* p = rho(i);
		const double* j = shapes.IPDets();
		const double* w = shapes.IPWeights();
		
		/* integrate */
		shapes.TopIP();
		while (shapes.NextIP()) {

			/* get shape function gradients */
			shapes.GradNa(grad_Na);
		
			/* integration factor */
			double jw = (*j++)*(*w++);
			double pm1_jw_by_V = ((*p++) - 1.0)*jw/V_0;
		
			/* loop over nodes */
			for (int k = 0; k < nodesX.Length(); k++) {
			
				int index = overlap_node_map.Map(nodesX[k]);
				if (index > -1) /* node is in overlap */ {
				
					/* inner product of bond and shape function gradient */
					double R_dot_dN = grad_Na.DotCol(k, R);
				
					/* assemble */
					f_a[index] += (R_dot_dN*pm1_jw_by_V);
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
	dArray2DT df_a_dp(f_a.Length(), df_dp.Length());
	df_a_dp = 0.0;
	dArray2DT df_dp_smoothing = df_dp;
	df_dp_smoothing = 0.0;
#if 0
	dArrayT diag(nip);
	diag = 0.0;
#endif

	/* compute residual */
	dArray2DT df(nen, nip);
	//dArrayT dp_i_dp(nip);
	dMatrixT ddp_i_dpdp(nip);
	df_dp = 0.0;
	ddf_dpdp = 0.0;
	dArrayT element_rho;
	dArrayT element_force;
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* set element information */
		const iArrayT& nodesX = coarse.ElementCard(overlap_cell[i]).NodesX();
		element_coords.SetLocal(nodesX);
		shapes.SetDerivatives();

		/* integration parameters */
		const double* p = rho(i);
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
			double pm1_jw = ((*p++) - 1.0)*jw;			
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
			df_dp(i,ip) += k2*pm1_jw;
			ddp_i_dpdp(ip,ip) += k2*jw;
					
			/* integrate the bond density term over the element */
			for (int k = 0; k < nodesX.Length(); k++) {
			
				int index = overlap_node_map.Map(nodesX[k]);
				if (index > -1) /* node is in overlap */ {
				
					/* inner product of bond and shape function gradient */
					double R_dot_dN = grad_Na.DotCol(k, R);
				
					/* assemble */
					df(k, ip) += (R_dot_dN*jw_by_V);
				}
			}
		}
		
		/* accumulate */
		for (int j = 0; j < nodesX.Length(); j++) {
			
			int index = overlap_node_map.Map(nodesX[j]);
			if (index > -1) /* node is in overlap */ {
			
				/* across all integration points */
				for (int k = 0; k < nip; k++) {
				
					/* force */
					df_dp(i,k) += f_a[index]*df(j,k);
					
					/* stiffness */
//					diag[k] += df(j,k)*df(j,k);
//					ddf_dpdp(i,k) += df(j,k)*df(j,k);
					
					//TEMP - accumulate
					df_a_dp(index, i*nip + k) += df(j,k);
				}
			}
		}
		
		/* regularization contribution to force and stiffness */
//		df_dp_smoothing.RowAlias(i, element_force);
		df_dp.RowAlias(i, element_force);
		rho.RowAlias(i, element_rho);		
		ATA_int.Multx(element_rho, element_force, 1.0, dMatrixT::kAccumulate);
//		for (int j = 0; j < nip; j++)
//			diag[j] += ATA_int(j,j);
//			ddf_dpdp(i,j) += ATA_int(j,j);

		/* penalty regularization */
		ATA_int += ddp_i_dpdp;
		
		ddf_dpdp.AddBlock(i*nip, i*nip, ATA_int);
	}

//cout << "df_dp: " << df_dp.no_wrap() << '\n';
//cout << "df_dp_smoothing: " << df_dp_smoothing.no_wrap() << '\n';
//df_dp += df_dp_smoothing;
//cout << "df_dp: " << df_dp.no_wrap() << '\n';

	//TEMP - add coupling term
	for (int i = 0; i < df_a_dp.MajorDim(); i++)
		ddf_dpdp.Outer(df_a_dp(i), df_a_dp(i), 1.0, dMatrixT::kAccumulate);

#if 0
	for (int i = 0; i < ddf_dpdp.Length(); i++)
		cout << ddf_dpdp[i] << '\n';
	cout.flush();
#endif
}

void FEManagerT_bridging::Compute_df_dp(const dArrayT& R, double V_0,
	const ArrayT<char>& cell_type, const InverseMapT& overlap_cell_map, const InverseMapT& overlap_node_map, 
	const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, LAdMatrixT& ddf_dpdp) const
{
	const char caller[] = "FEManagerT_bridging::Compute_df_dp";

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

//TEMP
fMainOut << "f_a =\n" << f_a << endl;
//TEMP

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

//TEMP
fMainOut << "df_a_dp =\n" << df_a_dp << endl;
//TEMP

	//TEMP - add coupling term
	for (int i = 0; i < df_a_dp.MajorDim(); i++)
		ddf_dpdp.Outer(df_a_dp(i), df_a_dp(i), 1.0, dMatrixT::kAccumulate);
}

void FEManagerT_bridging::Compute_df_dp(const dArrayT& R, double V_0, const ArrayT<char>& cell_type, 
	const InverseMapT& overlap_cell_map, const ArrayT<int>& overlap_node, const InverseMapT& overlap_node_map,
	const iArray2DT& cell_eq_i,
	const RaggedArray2DT<int>& inv_connects_i, const RaggedArray2DT<int>& inv_equations_i,
	const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, GlobalMatrixT& ddf_dpdp) const
{
	const char caller[] = "FEManagerT_bridging::Compute_df_dp";

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

//TEMP
fMainOut << "f_a =\n" << f_a << endl;
//TEMP

	/* gradient work space */
	const ParentDomainT& parent_domain = shapes.ParentDomain();	
	ArrayT<dMatrixT> ip_gradient(nip);
	dMatrixT jacobian_inv(nsd);
	dMatrixT A(nsd,nip), ATA(nip);
	ElementMatrixT ATA_int(nip, ElementMatrixT::kSymmetric);
	for (int i = 0; i < nip; i++) {
		ip_gradient[i].Dimension(nsd, nip);
		parent_domain.IPGradientTransform(i, ip_gradient[i]);
	}

	/* initialize return values */
	df_dp = 0.0;
	ddf_dpdp.Clear();

	/* regularization contributions to the force and stiffness matrix */
	dArray2DT df(nen, nip);
	dMatrixT ddp_i_dpdp(nip);
	dArrayT element_rho;
	dArrayT element_force;
	iArrayT eqnos;
	for (int i = 0; i < cell_type.Length(); i++)
		if (cell_type[i] == p_x) /* unknown bond density */ 
		{
			/* index within list of overlap cells */
			int overlap_cell_index = overlap_cell_map.Map(i);
			if (overlap_cell_index == -1) /* should all be within overlap */
				ExceptionT::GeneralFail(caller);
		
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
			shapes.TopIP();
			while (shapes.NextIP()) {

				/* integration factor */
				int ip = shapes.CurrIP();
				double jw = (*j++)*(*w++);
				double pm1_jw = ((*p++) - 1)*jw; /* penalization */
				double jw_by_V = jw/V_0;
	
				/* integrate density gradient matrix */
				parent_domain.DomainJacobian(element_coords, ip, jacobian_inv);
				jacobian_inv.Inverse();
				A.MultATB(jacobian_inv, ip_gradient[ip]);
				ATA.MultATB(A,A);
				ATA_int.AddScaled(smoothing*jw, ATA);
			
				/* add penalized force term */
				df_dp(overlap_cell_index,ip) += k2*pm1_jw;
				ddp_i_dpdp(ip,ip) += k2*jw;
			}

			/* regularization contribution to force */
			df_dp.RowAlias(overlap_cell_index, element_force);
			rho.RowAlias(overlap_cell_index, element_rho);		
			ATA_int.Multx(element_rho, element_force, 1.0, dMatrixT::kAccumulate);

			/* penalty regularization */
			ATA_int += ddp_i_dpdp;
		
			/* assemble stiffness contribution */
			cell_eq_i.RowAlias(overlap_cell_index, eqnos);
			ddf_dpdp.Assemble(ATA_int, eqnos);
		}

	/* work space */
	dArray2DT df_a_dp;
	nVariArray2DT<double> df_a_dp_man(0, df_a_dp, nip);
	ElementMatrixT df_a_dp_2(ElementMatrixT::kSymmetric);
	nVariMatrixT<double> df_a_dp_2_man(0, df_a_dp_2);

	/* add Cauchy-Born contribution from "connectivities" of overlap nodes */
	iArrayT node_elements;
	for (int i = 0; i < inv_connects_i.MajorDim(); i++) {
	
		int node = overlap_node[i];
		int node_index = overlap_node_map.Map(node);
		if (node_index == -1) ExceptionT::GeneralFail(caller);
	
		/* elements in the nodal support */
		inv_connects_i.RowAlias(i, node_elements);
	
		/* dimension workspace */
		df_a_dp_man.SetMajorDimension(inv_connects_i.MinorDim(i), false);
		
		/* initialize */
		df_a_dp = 0.0;
		for (int e = 0; e < node_elements.Length(); e++) {

			int element = node_elements[e];

			/* index within list of overlap cells */
			int overlap_cell_index = overlap_cell_map.Map(element);
			if (overlap_cell_index == -1) /* should all be within overlap */
				ExceptionT::GeneralFail(caller);
		
			/* set element information */
			const iArrayT& nodesX = coarse->ElementCard(element).NodesX();
			element_coords.SetLocal(nodesX);
			shapes.SetDerivatives();
			
			/* find the local number of the node within the element */
			int local_node = -1;
			for (int k = 0; local_node == -1 && k < nodesX.Length(); k++)
				if (nodesX[k] == node)
					local_node = node;
			if (local_node == -1) ExceptionT::GeneralFail(caller);

			/* integration parameters */
			const double* p = rho(overlap_cell_index);
			const double* j = shapes.IPDets();
			const double* w = shapes.IPWeights();

			shapes.TopIP();
			while (shapes.NextIP()) {

				/* integration factor */
				int ip = shapes.CurrIP();
				double jw = (*j++)*(*w++);
				double jw_by_V = jw/V_0;
	
				/* get shape function gradients */
				shapes.GradNa(grad_Na);

				/* inner product of bond and shape function gradient */
				double R_dot_dN = grad_Na.DotCol(local_node, R);
				
				/* assemble */
				df_a_dp(e, ip) += R_dot_dN*jw_by_V;
			}		
		}

//TEMP
fMainOut << "df_a_dp =\n" << df_a_dp << endl;
//TEMP

		/* assemble stiffness */
		df_a_dp_2_man.SetDimensions(inv_equations_i.MinorDim(i));
		df_a_dp_2.Outer(df_a_dp, df_a_dp);
		inv_equations_i.RowAlias(i, eqnos);
		ddf_dpdp.Assemble(df_a_dp_2, eqnos);
		
		/* force contribution */
		for (int j = 0; j < eqnos.Length(); j++)
			df_dp[eqnos[j]-1] += f_a[node_index]*df_a_dp[j];
	}
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
	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail("FEManagerT_bridging::GhostNodeBonds", "interpolation data not set");
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
}

void FEManagerT_bridging::GhostNodeBonds(const dArrayT& R_i, const dArray2DT& point_coords, 
	const RaggedArray2DT<int>& ghost_neighbors_all, RaggedArray2DT<int>& ghost_neighbors_i, 
	AutoArrayT<int>& overlap_cell_i, AutoArrayT<int>& overlap_node_i) const
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

#ifdef __DEBUG__
	AutoArrayT<int> a0;
	AutoArrayT<int> a1;
#endif

	/* keep only bonds (from free atoms) to ghost points */
	int num_overlap = 0;
	iArrayT neighbors_i;
	AutoFill2DT<int> gh_neigh(ghost_neighbors_all.MajorDim(), 1, 25, 10);
	dArrayT bond(R_i.Length());
	double R = R_i.Magnitude();
	for (int i = 0; i < ghost_neighbors_all.MajorDim(); i++) {

		/* the neighbor list */
		ghost_neighbors_all.RowAlias(i, neighbors_i);
		
		/* self is first neighbor */
		gh_neigh.Append(i, neighbors_i[0]);
			
		/* search other neighbors for ghosts */
		for (int j = 1; j < neighbors_i.Length(); j++) {
		
			/* bond terminating at follower node */
			bond.DiffOf(point_coords(neighbors_i[j]), point_coords(neighbors_i[0]));
			double L_b = bond.Magnitude();
			
			/* bond direction and length */
			double cosRR = dArrayT::Dot(R_i, bond)/L_b/R;
			if (fabs(R - L_b)/R < 1.0e-02 && fabs(fabs(cosRR) - 1.0) < kSmall) {
			
				/* keep neighbor */
				gh_neigh.Append(i, neighbors_i[j]);

#ifdef __DEBUG__
				/* collect bond atoms */
				a0.Append(neighbors_i[0]);
				a1.Append(neighbors_i[j]);
#endif

				/* mark cell */
				int neighbor_local = follower_point_map.Map(neighbors_i[j]);	
				int cell = interpolating_cell[neighbor_local];
				if (is_overlap_cell[cell] == 'f') {
					num_overlap++;
					is_overlap_cell[cell] = 't';
				}
			}
		}
	}

	/* copy/compress */
	ghost_neighbors_i.Copy(gh_neigh);

#ifdef __DEBUG__
	fMainOut << "\n contributing bonds:\n";
	for (int i = 0; i < a0.Length(); i++)
		fMainOut << a0[i]+1 << " " << a1[i]+1 << '\n';
	fMainOut.flush();
#endif
	
	/* collect overlap cells */
	overlap_cell_i.Dimension(num_overlap);
	num_overlap = 0;
	for (int i = 0; i < is_overlap_cell.Length(); i++)
		if (is_overlap_cell[i] == 't')
			overlap_cell_i[num_overlap++] = i;

	/* mark overlap nodes */
	int nnd = fNodeManager->NumNodes();	
	ArrayT<char> is_overlap_node(nnd);
	is_overlap_node = 'f';
	num_overlap = 0;
	for (int i = 0; i < overlap_cell_i.Length(); i++) {
		const iArrayT& nodesX = coarse->ElementCard(overlap_cell_i[i]).NodesX();
		for (int j = 0; j < nodesX.Length(); j++) {
			int nd = nodesX[j];
			if (is_overlap_node[nd] == 'f') {
				num_overlap++;
				is_overlap_node[nd] = 't';
			}
		}
	}

	/* collect overlap nodes */
	overlap_node_i.Dimension(num_overlap);
	num_overlap = 0;
	for (int i = 0; i < is_overlap_node.Length(); i++)
		if (is_overlap_node[i] == 't')
			overlap_node_i[num_overlap++] = i;
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
		int active_node_index = active_node_map.Map(nodes[i]);
		if (active_node_index != -1)
			elements_per_node.Append(active_node_index, active_elements[i]);
	}
	
	/* copy into return value */
	transpose_connects.Copy(elements_per_node);
}

#endif  /* BRIDGING_ELEMENT */
