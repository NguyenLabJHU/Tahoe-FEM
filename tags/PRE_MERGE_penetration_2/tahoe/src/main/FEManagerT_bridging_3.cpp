/* $Id: FEManagerT_bridging_3.cpp,v 1.2 2004-06-26 06:24:58 paklein Exp $ */
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
#include "EAMFCC3D.h"
#include "EAMT.h"
#include "ShapeFunctionT.h"
#include "dSymMatrixT.h"
#include "ElementSupportT.h"

/* headers needed to compute the correction for overlap */
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

void FEManagerT_bridging::CorrectOverlap_3(const RaggedArray2DT<int>& point_neighbors, const dArray2DT& point_coords, 
	double smoothing, double k2, int nip)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_3";

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
	if (nip != 1 && nip != coarse->NumIP())
		ExceptionT::GeneralFail(caller, "number of integration points needs to be 1 or %d: %d",
			coarse->NumIP(), nip);

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
	dArray2DT bond_densities(overlap_cell_all.Length(), nip*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space that changes for each bond family */
	CCSMatrixT ddf_dpdp_i(Output(), GlobalMatrixT::kZeroPivots);
	dArrayT rhs;
	VariArrayT<double> rhs_man(0, rhs);
	
	dArray2DT p_i, dp_i, df_dp_i;
	nArray2DGroupT<double> ip_unknown_group(0, false, nip);
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
	iArray2DT bond_densities_i_eq_all;
	nVariArray2DT<int> bond_densities_i_eq_all_man(0, bond_densities_i_eq_all, nip);
	iArray2DT bond_densities_i_eq_active;
	nVariArray2DT<int> bond_densities_i_eq_active_man(0, bond_densities_i_eq_active, nip);
	RaggedArray2DT<int> inv_connects_i; 
	RaggedArray2DT<int> inv_equations_all_i;
	RaggedArray2DT<int> inv_equations_active_i;
	AutoArrayT<int> bondfree_cell_i;
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond family at a time */ {

		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
	
		/* collect "acive" ghost node bonds - bonds of type R_i terminating at a ghost node */
		RaggedArray2DT<int> ghost_neighbors_i;
		GhostNodeBonds_2(R_i, point_coords, ghost_neighbors_all, ghost_neighbors_i, overlap_cell_i, overlap_node_i);
		overlap_cell_i_map.SetMap(overlap_cell_i);
		overlap_node_i_map.SetMap(overlap_node_i);
		if (fPrintInput) {
			iArrayT tmp;
			
			fMainOut << "overlap nodes for bond: " << i+1 << ": {" << R_i.no_wrap() << "}:\n";
			tmp.Alias(overlap_node_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;

			fMainOut << "overlap cells for bond: " << i+1 << ": {" << R_i.no_wrap() << "}:\n";
			tmp.Alias(overlap_cell_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;
		}
		
		/* collect list of cells not containing any active bonds */
		BondFreeElements(ghost_neighbors_i, bondfree_cell_i);

		/* classify element types */
		cell_type = p_0; /* all zero density */
		for (int j = 0; j < bondfree_cell_i.Length(); j++) /* full density */
			cell_type[bondfree_cell_i[j]] = p_1;
		for (int j = 0; j < overlap_cell_i.Length(); j++) /* unknown density */
			cell_type[overlap_cell_i[j]] = p_x;
			
		/* dimension work space */
		ip_unknown_group.SetMajorDimension(overlap_cell_i.Length(), false);
		overlap_node_group.Dimension(overlap_node_i.Length(), false);
		
		/* number unknown bond densities */
		bond_densities_i_eq_all_man.SetMajorDimension(overlap_cell_i.Length(), false);
		bond_densities_i_eq_active_man.SetMajorDimension(overlap_cell_i.Length(), false);
		CountIPBonds(overlap_cell_i, bond_densities_i_eq_active);
		int num_eq_active = 0;
		int num_eq_all = 0;
		for (int j = 0; j < overlap_cell_i.Length(); j++)
			for (int k = 0; k < bond_densities_i_eq_all.MinorDim(); k++) {

				/* mark all densities */
				bond_densities_i_eq_all(j,k) = ++num_eq_all;				

				/* mark unknown densities */
				if (bond_densities_i_eq_active(j,k) > 0)
					bond_densities_i_eq_active(j,k) = ++num_eq_active;
				else
					bond_densities_i_eq_active(j,k) = -1; /* inactive */
		}
		rhs_man.SetLength(num_eq_active, false);

		/* "inverse" connectivities for overlap cells in the support of overlap nodes */
		TransposeConnects(*coarse, overlap_node_i, overlap_cell_i, inv_connects_i);
		
		/* collect equation numbers for "inverse" elements */
		inv_equations_all_i.Configure(inv_connects_i, nip);
		inv_equations_active_i.Configure(inv_connects_i, nip);
		for (int j = 0; j < inv_connects_i.MajorDim(); j++) {
			const int* element_per_node = inv_connects_i(j);
			int* equations_per_node_all = inv_equations_all_i(j);
			int* equations_per_node_active = inv_equations_active_i(j);
			for (int k = 0; k < inv_connects_i.MinorDim(j); k++) {
				int overlap_cell_index = overlap_cell_i_map.Map(element_per_node[k]);
				int* equations_all = bond_densities_i_eq_all(overlap_cell_index);
				int* equations_active = bond_densities_i_eq_active(overlap_cell_index);
				for (int l = 0; l < bond_densities_i_eq_all.MinorDim(); l++) {
					*equations_per_node_all++ = *equations_all++;
					*equations_per_node_active++ = *equations_active++;
				}
			}
		}

		/* configure linear solver */
		ddf_dpdp_i.AddEquationSet(inv_equations_active_i);
		ddf_dpdp_i.Initialize(num_eq_active, num_eq_active, 1);
		
		/* compute contribution from bonds terminating at "ghost" atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors_i, point_coords, overlap_node_i_map, sum_R_N);
		if (fPrintInput) {
			fMainOut << "ghost bond contritbution:\n";
			fMainOut << "R.sum_R_N =\n" << sum_R_N << endl;
		}

		/* initialize */
		p_i = 1.0;
		
		/* compute residual - add Cauchy-Born contribution */
		f_a = sum_R_N;
		Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
			bond_densities_i_eq_active, inv_connects_i, inv_equations_all_i, inv_equations_active_i,
			p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
		if (fPrintInput) {
			fMainOut << "residual =\n" << df_dp_i << endl;
		}

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

				/* copy to rhs */
				for (int j = 0; j < bond_densities_i_eq_active.Length(); j++) {
					int eq = bond_densities_i_eq_active[j] - 1;
					if (eq >= 0) /* active */
						rhs[eq] = -df_dp_i[j];
				}

				/* solve system */
				ddf_dpdp_i.Solve(rhs);

#if 0
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
#endif
				/* update densities */
				for (int j = 0; j < bond_densities_i_eq_active.Length(); j++) {
					int eq = bond_densities_i_eq_active[j] - 1;
					if (eq >= 0) /* active */
						p_i[j] += rhs[eq];
				}

				/* recompute residual */			
				f_a = sum_R_N;
				Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map,
					bond_densities_i_eq_active, inv_connects_i, inv_equations_all_i, inv_equations_active_i,
					p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
				
//				Compute_df_dp(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
				error = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
				cout << "iteration = " << iter 
				     << "\n e/e_0 = " << error/error_0 
				     << "\n||fa|| = " << sqrt(dArrayT::Dot(f_a, f_a)) << endl;

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
			for (int j = 0; j < nip; j++)
				bond_densities(cell, i+j*nb) = p_i(k,j);
		}		
	}
	
	/* write densities */
	if (fPrintInput) {
	
		/* dimensions */
		const NodeManagerT* node_manager = NodeManager();
		int nsd = node_manager->NumSD();
		int nen = coarse->NumElementNodes();

		/* shape functions */
		LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
		NodeManager()->RegisterCoordinates(element_coords);
		ShapeFunctionT shapes = ShapeFunctionT(coarse->GeometryCode(), nip, element_coords);
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
	
	/* bound bond densities [0,1] */
	for (int i = 0; i < bond_densities.Length(); i++)
		if (bond_densities[i] > 1.0)
			bond_densities[i] = 1.0;
		else if (bond_densities[i] < 0.0)
			bond_densities[i] = 0.0;

	/* write unknowns into the state variable space */
	ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
	if (nip == coarse->NumIP()) {
		for (int i = 0; i < overlap_cell_all.Length(); i++) {
	
			/* element information */
			ElementCardT& element = non_const_coarse->ElementCard(overlap_cell_all[i]);
	
			/* allocate space */
			element.Dimension(0, bond_densities.MinorDim());
		
			/* copy in densities */
			element.DoubleData() = bond_densities(i);
		}
	}
	else if (nip == 1) /* constant across all integration points */ {
		int num_elem_ip = coarse->NumIP();
		dArray2DT ip_bond_densities;
		for (int i = 0; i < overlap_cell_all.Length(); i++) {
	
			/* element information */
			ElementCardT& element = non_const_coarse->ElementCard(overlap_cell_all[i]);
	
			/* allocate space */
			element.Dimension(0, num_elem_ip*bond_densities.MinorDim());
		
			/* copy in densities */
			ip_bond_densities.Alias(num_elem_ip, bond_densities.MinorDim(), element.DoubleData().Pointer());
			for (int j = 0; j < num_elem_ip; j++)
				ip_bond_densities.SetRow(j, bond_densities(i));
		}	
	}
	else ExceptionT::GeneralFail(caller);	
}

#endif  /* BRIDGING_ELEMENT */
