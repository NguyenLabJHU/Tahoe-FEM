/* $Id: FEManagerT_bridging_4.cpp,v 1.2 2004-06-26 06:24:58 paklein Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "CommManagerT.h"

#include "BridgingScaleT.h"
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
#include "SecantMethodT.h"

#include "CCSMatrixT.h"
#ifdef __SPOOLES__
#include "SPOOLESMatrixT.h"
#endif

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

void FEManagerT_bridging::CorrectOverlap_4(const RaggedArray2DT<int>& point_neighbors, const dArray2DT& point_coords, 
	double smoothing, double k2, double k_r, int nip)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_4";

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
	dArray2DT bonds = bond_lattice.Bonds();
	double V_0 = (hex_2D) ? hex_2D->CellVolume() : ((fcc_3D) ? fcc_3D->CellVolume() : mat_1D->CellVolume());
	double R_0 = (hex_2D) ? hex_2D->NearestNeighbor() : ((fcc_3D) ? fcc_3D->NearestNeighbor() : mat_1D->NearestNeighbor());
	bonds *= R_0;

	/* group bonds in shells */
	iArrayT shell; /* shell number of each bond */
	dArrayT shell_bond_length; /* length of bonds in each shell */
	shell = -1;
	int num_shells = NumberShells(bonds, shell, shell_bond_length);

	/* unknown bond densities */
	dArray2DT bond_densities(overlap_cell_all.Length(), nip*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space that changes for each shell */
#ifdef __SPOOLES__
	SPOOLESMatrixT ddf_dpdp_i(Output(), GlobalMatrixT::kZeroPivots, true, true);
#else
	CCSMatrixT ddf_dpdp_i(Output(), GlobalMatrixT::kZeroPivots);
#endif	
	
	dArray2DT p_i, dp_i, df_dp_i;
	nArray2DGroupT<double> ip_unknown_group(0, false, nip);
	ip_unknown_group.Register(p_i);
	ip_unknown_group.Register(dp_i);
	ip_unknown_group.Register(df_dp_i);
	dArray2DT sum_R_N;
	dArray2DT f_a;
	nArray2DGroupT<double> overlap_node_group(0, false, point_coords.MinorDim());
	overlap_node_group.Register(sum_R_N);
	overlap_node_group.Register(f_a);

	/* constraints */
	iArrayT eqnos;
	iArray2DT constraint_eq;
	nVariArray2DT<int> constraint_eq_man(0, constraint_eq, 2);
	dArrayT f_constraint(2);
	ElementMatrixT K_constraint(2, ElementMatrixT::kSymmetric);
	dArrayT rhs;
	VariArrayT<double> rhs_man(0, rhs);
	dArrayT update;
	VariArrayT<double> update_man(0, update);
	dArrayT constraint;
	VariArrayT<double> constraint_man(0, constraint);

	/* secant method solver */
	SecantMethodT secant_search(10, 0.1);

	/* solve unknowns */
	int bond_density_offset = 0;
	ArrayT<char> cell_type(nel);
	AutoArrayT<int> overlap_cell_i;
	InverseMapT overlap_cell_i_map;
	overlap_cell_i_map.SetOutOfRange(InverseMapT::MinusOne);
	AutoArrayT<int> overlap_node_i;
	InverseMapT overlap_node_i_map;
	overlap_node_i_map.SetOutOfRange(InverseMapT::MinusOne);
	iArray2DT bond_densities_i_eq;
	nVariArray2DT<int> bond_densities_i_eq_man(0, bond_densities_i_eq, nip);
	RaggedArray2DT<int> inv_connects_i; 
	RaggedArray2DT<int> inv_equations_i; 
	AutoArrayT<int> bondfree_cell_i;
	for (int i = 0; i < shell_bond_length.Length(); i++) /* one shell at a time */ {

		/* count number of bonds in the shell */
		int num_shell_bonds = 0;
		for (int j = 0; j < shell.Length(); j++)
			if (shell[j] == i)
				num_shell_bonds++;

		/* collect shell bonds */
		dArray2DT shell_bonds(num_shell_bonds, bonds.MinorDim());
		num_shell_bonds = 0;
		for (int j = 0; j < shell.Length(); j++)
			if (shell[j] == i)
				shell_bonds.SetRow(num_shell_bonds++, bonds(j));

		/* collect "acive" ghost node bonds - bonds in same shell terminating at a ghost node */
		RaggedArray2DT<int> ghost_neighbors_i;
		GhostNodeBonds_4(shell_bond_length[i], point_coords, ghost_neighbors_all, ghost_neighbors_i, overlap_cell_i, overlap_node_i);
		overlap_cell_i_map.SetMap(overlap_cell_i);
		overlap_node_i_map.SetMap(overlap_node_i);
		if (fPrintInput) {
			iArrayT tmp;
			
			fMainOut << "overlap nodes for shell: " << i+1 << ": " << shell_bond_length[i] << ":\n";
			tmp.Alias(overlap_node_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;

			fMainOut << "overlap cells for shell: " << i+1 << ": " << shell_bond_length[i] << ":\n";
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
		ip_unknown_group.Dimension(overlap_cell_i.Length(), num_shell_bonds*nip);		
		overlap_node_group.SetMajorDimension(overlap_node_i.Length(), false);
		
		/* number unknown bond densities */
		bond_densities_i_eq_man.Dimension(overlap_cell_i.Length(), num_shell_bonds*nip);
		bond_densities_i_eq = -1;
		int num_eq = 0;
		for (int j = 0; j < bond_densities_i_eq.MajorDim(); j++)
			for (int k = 0; k < bond_densities_i_eq.MinorDim(); k++)
				bond_densities_i_eq(j,k) = ++num_eq;
		int num_eq_p = num_eq;

		/* "inverse" connectivities for overlap cells in the support of overlap nodes */
		TransposeConnects(*coarse, overlap_node_i, overlap_cell_i, inv_connects_i);
		
		/* collect equation numbers for "inverse" elements */
		inv_equations_i.Configure(inv_connects_i, nip*num_shell_bonds);
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
		
		/* constraint equations - one constraint per unknown */
		int num_eq_L = 0;
		if (k_r > kSmall) {
		
			/* dimension equations array */
			num_eq_L = num_eq_p;
			constraint_eq_man.SetMajorDimension(num_eq_L, false);
		
			/* assign equation numbers */
			for (int j = 0; j < bond_densities_i_eq.Length(); j++) {
				constraint_eq(j,0) = bond_densities_i_eq[j];
				constraint_eq(j,1) = ++num_eq;
			}
		
			/* configure solver */
			ddf_dpdp_i.AddEquationSet(constraint_eq);

			/* initialize unknowns */
			constraint_man.SetLength(num_eq_L, false);
			constraint = 0.0;
		}

		/* configure equation system */
		ddf_dpdp_i.AddEquationSet(inv_equations_i);
		ddf_dpdp_i.Initialize(num_eq, num_eq, 1);
		rhs_man.SetLength(num_eq, false);
		update_man.SetLength(num_eq, false);		

		/* compute contribution from bonds terminating at "ghost" atoms */
		ComputeSum_signR_Na_4(shell_bond_length[i], ghost_neighbors_i, point_coords, overlap_node_i_map, sum_R_N);
		if (fPrintInput) {
			fMainOut << "ghost bond contritbution:\n";
			fMainOut << "sum_R_N =\n" << sum_R_N << endl;
		}

		/* initialize */
		p_i = 1.0;
		rhs = 0.0;
		
		/* compute residual - add Cauchy-Born contribution */
		f_a = sum_R_N;
		Compute_df_dp_4(shell_bonds, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
			bond_densities_i_eq, inv_connects_i, inv_equations_i,
			p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
			
		/* assemble residual */
		for (int j = 0; j < df_dp_i.Length(); j++)
			rhs[j] -= df_dp_i[j];
			
		/* apply constraints */
		for (int j = 0; j < num_eq_L; j++)
		{
			/* equation numbers */
			constraint_eq.RowAlias(j, eqnos);

			/* constraint equation: g >= 0 */
			double   g = 0.25 - (p_i[j] - 0.5)*(p_i[j] - 0.5);
			double  Dg = -2.0*(p_i[j] - 0.5);
			double DDg = -2.0;
		
			/* augmented multiplier */
			double L_r = constraint[j] + k_r*g;
			if (L_r < 0.0) /* active constraint */
			{
				/* force */
				rhs[eqnos[0]-1] -= Dg*L_r;
				rhs[eqnos[1]-1] -= g;

				/* stiffness */
				K_constraint(0,0) = DDg*L_r + k_r*Dg*Dg;
				K_constraint(0,1) = Dg;
				K_constraint(1,0) = Dg;
				K_constraint(1,1) = 0.0;
			}
			else /* inactive */
			{
				/* force */
				rhs[eqnos[0]-1] += 0.0;
				rhs[eqnos[1]-1] += constraint[j]/k_r;
			
				/* stiffness */
				K_constraint(0,0) = 0.0;
				K_constraint(0,1) = 0.0;
				K_constraint(1,0) = 0.0;
				K_constraint(1,1) = -1.0/k_r;
			}
			
			/* assembly stiffness constribution */
			ddf_dpdp_i.Assemble(K_constraint, eqnos);
		}
			
		if (fPrintInput) {
			fMainOut << "residual =\n" << df_dp_i << endl;
		}

		/* solve bond densities */
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+06;		
		int max_iter = 200;
		int iter = 0;
		double error_0 = sqrt(dArrayT::Dot(rhs, rhs));		
		double error = error_0;
		cout << "\n{e, ||fa||} = {" << error << ", " <<  sqrt(dArrayT::Dot(f_a, f_a)) << "}" << endl;
		while (iter++ < max_iter && error > abs_tol && error/error_0 > rel_tol && error/error_0 < div_tol) {

			/* catch errors in linear solver */
			try {

				/* solve system */
				update = rhs;
				ddf_dpdp_i.Solve(update);

				/* line search */
				int ls_count = 0;
				bool ls_continue = true;
				double s = 0.0;
				double s_last = 0.0;
				double G_R_0 = 0.0;
				while (ls_continue)
				{
					/* select step size */
					if (ls_count == 0)
						s = 1.0; /* full Newton step */
					else if (ls_count == 1)
						s = 0.5;
					else /* secant method */
						s = secant_search.NextGuess();

					/* scale update */
					rhs = update;
					rhs *= (s - s_last);

					/* update densities */
					for (int j = 0; j < num_eq_p; j++)
						p_i[j] += rhs[j];
								
					/* update constraints */
					for (int j = 0; j < num_eq_L; j++)
						constraint[j] += rhs[j+num_eq_p];

					/* recompute residual */			
					f_a = sum_R_N;
					Compute_df_dp_4(shell_bonds, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map,
						bond_densities_i_eq, inv_connects_i, inv_equations_i,
						p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
					
					/* assemble residual */
					rhs = 0.0;
					for (int j = 0; j < df_dp_i.Length(); j++)
						rhs[j] -= df_dp_i[j];
				
					/* apply constraints */
					for (int j = 0; j < num_eq_L; j++)
					{
						/* equation numbers */
						constraint_eq.RowAlias(j, eqnos);

						/* constraint equation: g >= 0 */
						double   g = 0.25 - (p_i[j] - 0.5)*(p_i[j] - 0.5);
						double  Dg = -2.0*(p_i[j] - 0.5);
						double DDg = -2.0;
		
						/* augmented multiplier */
						double L_r = constraint[j] + k_r*g;
						if (L_r < 0.0) /* active constraint */
						{
							/* force */
							rhs[eqnos[0]-1] -= Dg*L_r;
							rhs[eqnos[1]-1] -= g;

							/* stiffness */
							K_constraint(0,0) = DDg*L_r + k_r*Dg*Dg;
							K_constraint(0,1) = Dg;
							K_constraint(1,0) = Dg;
							K_constraint(1,1) = 0.0;
						}
						else /* inactive */
						{
							/* force */
							rhs[eqnos[0]-1] += 0.0;
							rhs[eqnos[1]-1] += constraint[j]/k_r;
			
							/* stiffness */
							K_constraint(0,0) = 0.0;
							K_constraint(0,1) = 0.0;
							K_constraint(1,0) = 0.0;
							K_constraint(1,1) = -1.0/k_r;
						}
			
						/* assembly stiffness constribution */
						ddf_dpdp_i.Assemble(K_constraint, eqnos);
					}
					
					/* orthogonality */
					double G_R = dArrayT::Dot(rhs, update);

					/* process results */
					if (ls_count == 0)
					{
						G_R_0 = G_R;
					
						/* quick exit */
						if (fabs(G_R_0) < kSmall)
							ls_continue = false;
					}
					else if (ls_count == 1) /* initialize secant search */
						secant_search.Reset(s, G_R, 1.0, G_R_0);
					else
					{
						int test = secant_search.NextPoint(s, G_R);
						if (test == 1)
							ls_continue = false;
						else if (test == -1)
							ExceptionT::GeneralFail(caller, "line search failed");
					}

					/* next iteration */
					ls_count++;
					s_last = s;
				}
					
				error = sqrt(dArrayT::Dot(rhs, rhs));
				cout << setw(5) << iter << ": {n_ls, s, e/e_0, ||fa||} = {"<< secant_search.Iterations() << ", " << s << ", " << error/error_0 << ", " <<  sqrt(dArrayT::Dot(f_a, f_a)) << "}" << endl;
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
			for (int b = 0; b < num_shell_bonds; b++)
				for (int j = 0; j < nip; j++)
					bond_densities(cell, (bond_density_offset+b) + j*nb) = p_i(k, b*nip + j);
		}

		/* next shell of bonds */
		bond_density_offset += num_shell_bonds;
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
	
#if 0
	/* bound bond densities [0,1] */
	for (int i = 0; i < bond_densities.Length(); i++)
		if (bond_densities[i] > 1.0)
			bond_densities[i] = 1.0;
		else if (bond_densities[i] < 0.0)
			bond_densities[i] = 0.0;
#endif

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

void FEManagerT_bridging::GhostNodeBonds_4(double R_i, const dArray2DT& point_coords, 
	const RaggedArray2DT<int>& ghost_neighbors_all, RaggedArray2DT<int>& ghost_neighbors_i, 
	AutoArrayT<int>& overlap_cell_i, AutoArrayT<int>& overlap_node_i) const
{
	const char caller[] = "FEManagerT_bridging::GhostNodeBonds_4";

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
	dArrayT bond(point_coords.MinorDim());
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
			if (fabs(R_i - L_b)/R_i < 1.0e-02) {
			
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

	/* mark overlap nodes - nodes whose support contains an active
	 * ghost atom bond */
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

/* compute contribution from bonds to ghost atoms */
void FEManagerT_bridging::ComputeSum_signR_Na_4(double R_i, const RaggedArray2DT<int>& ghost_neighbors, 
	const dArray2DT& coords, const InverseMapT& overlap_node_map, dArray2DT& sum_R_N) const
{
	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail("FEManagerT_bridging::ComputeSum_signR_Na_4", "interpolation data not set");

	/* interpolating data */
	const RaggedArray2DT<int>& point_in_cell = fFollowerCellData.PointInCell();
	const InverseMapT& follower_point_map = fFollowerCellData.GlobalToLocal();
	const iArrayT& interpolating_cell = fFollowerCellData.InterpolatingCell();
	const dArray2DT& interpolating_weights = fFollowerCellData.InterpolationWeights();

	/* accumulate ghost bond contribution */
	sum_R_N = 0.0;
	dArrayT bond(coords.MinorDim());
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
				if (fabs(R_i - L_b)/R_i < 1.0e-02) {

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
							sum_R_N.AddToRowScaled(overlap_node_index, weights[k], bond);
					}
				}		
			}
		}
}

void FEManagerT_bridging::Compute_df_dp_4(const dArray2DT& shell_bonds, double V_0, const ArrayT<char>& cell_type, 
	const InverseMapT& overlap_cell_map, const ArrayT<int>& overlap_node, const InverseMapT& overlap_node_map,
	const iArray2DT& cell_eq_i, const RaggedArray2DT<int>& inv_connects_i, const RaggedArray2DT<int>& inv_equations_i,
	const dArray2DT& rho, dArray2DT& f_a, double smoothing, double k2, dArray2DT& df_dp, GlobalMatrixT& ddf_dpdp) const
{
	const char caller[] = "FEManagerT_bridging::Compute_df_dp_4";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");

	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse->NumElementNodes();
	int nsh = shell_bonds.MajorDim();
	int nip = rho.MinorDim()/nsh;

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse->GeometryCode(), nip, element_coords);
	shapes.Initialize();

	/* integrate bond density term */
	dArrayT R_i;
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

			/* loop over shell bonds */
			for (int b = 0; b < nsh; b++)
			{
				/* bond */
				shell_bonds.RowAlias(b, R_i);
			
				/* bond density */
				const double* p = NULL;
				if (cell_type[i] == p_x) {
					int overlap_cell_index = overlap_cell_map.Map(i);
					if (overlap_cell_index == -1) ExceptionT::GeneralFail(caller);
					p = rho(overlap_cell_index) + b*nip;
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
							double R_dot_dN = grad_Na.DotCol(k, R_i);
				
							/* assemble */
							f_a.AddToRowScaled(overlap_node_index, R_dot_dN*p_jw_by_V, R_i);
						}
					}
				}
			}			
		}
	if (fPrintInput) {
		fMainOut << "f_a =\n" << f_a << endl;
	}

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

			/* loop over bonds in the shell */
			for (int b = 0; b < nsh; b++)
			{
				/* integration parameters */
				const double* p = rho(overlap_cell_index) + b*nip;
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
					df_dp(overlap_cell_index, b*nip + ip) += k2*pm1_jw;
					ddp_i_dpdp(ip,ip) += k2*jw;
				}

				/* regularization contribution to force */
				element_force.Alias(nip, df_dp(overlap_cell_index) + b*nip);
				element_rho.Alias(nip, rho(overlap_cell_index) + b*nip);
				ATA_int.Multx(element_rho, element_force, 1.0, dMatrixT::kAccumulate);

				/* penalty regularization */
				ATA_int += ddp_i_dpdp;
		
				/* assemble stiffness contribution */
				eqnos.Alias(nip, cell_eq_i(overlap_cell_index) + b*nip);
				ddf_dpdp.Assemble(ATA_int, eqnos);
			}
		}

	/* work space */
	dMatrixT df_a_dp;
	nVariMatrixT<double> df_a_dp_man(0, df_a_dp);
	dArrayT dp(nsh*nip);
	VariArrayT<double> dp_man(0, dp);
	ElementMatrixT df_a_dp_2(ElementMatrixT::kSymmetric);
	nVariMatrixT<double> df_a_dp_2_man(0, df_a_dp_2);

	/* add Cauchy-Born contribution from "connectivities" of overlap nodes */
	iArrayT node_elements;
	dArrayT f_a_i;
	for (int i = 0; i < inv_connects_i.MajorDim(); i++) {
	
		int node = overlap_node[i];
		int node_index = overlap_node_map.Map(node);
		if (node_index == -1) ExceptionT::GeneralFail(caller);
	
		/* elements in the nodal support */
		inv_connects_i.RowAlias(i, node_elements);
	
		/* dimension workspace */
		df_a_dp_man.SetDimensions(nsd, node_elements.Length()*nsh*nip);
		dp_man.SetLength(node_elements.Length()*nsh*nip, false);
		
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
					local_node = k;
			if (local_node == -1) ExceptionT::GeneralFail(caller);

			/* loop over shell bonds */
			for (int b = 0; b < nsh; b++)
			{
				/* bond */
				shell_bonds.RowAlias(b, R_i);
					
				/* integration parameters */
				const double* p = rho(overlap_cell_index) + b*nip;
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
					double R_dot_dN = grad_Na.DotCol(local_node, R_i);

					/* assemble */
					for (int l = 0; l < nsd; l++)
						df_a_dp(l, e*nsh*nip + b*nip + ip) += R_i[l]*R_dot_dN*jw_by_V;
				}
			}		
		}

		/* assemble stiffness */
		df_a_dp_2_man.SetDimensions(inv_equations_i.MinorDim(i));
		df_a_dp_2.MultATB(df_a_dp, df_a_dp);
		inv_equations_i.RowAlias(i, eqnos);
		ddf_dpdp.Assemble(df_a_dp_2, eqnos);
		
		/* force contribution */
		f_a.RowAlias(node_index, f_a_i);
		df_a_dp.MultTx(f_a_i, dp);
		inv_equations_i.RowAlias(i, eqnos);
		for (int j = 0; j < eqnos.Length(); j++) {
			df_dp[eqnos[j]-1] += dp[j];
		}
	}
}

#endif  /* BRIDGING_ELEMENT */
