/* $Id: FEManagerT_bridging_22.cpp,v 1.10 2005-02-18 02:37:21 paklein Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "CommManagerT.h"

#include "ifstreamT.h"
#include "BridgingScaleT.h"
#include "ParticleT.h"
#include "ParticlePairT.h"
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
#include "SecantMethodT.h"

#include "CCSMatrixT.h"
#ifdef __SPOOLES__
#include "SPOOLESMatrixT.h"
#endif

#include <float.h>

/* debugging */
//#define __DEBUG__ 1
#undef __DEBUG__

/* atom/point types */
const char free_ = 'f';
const char not_free_ = 'n';

/* element types */
const char p_0 = 'a'; /* bond density = 0 */
const char p_1 = 'b'; /* bond density = 1 */
const char p_x = 'c'; /* unknown: 0 < bond density < 1 */

const double sqrt2 = sqrt(2.0);

using namespace Tahoe;

void FEManagerT_bridging::CorrectOverlap_22(const RaggedArray2DT<int>& point_neighbors, const dArray2DT& point_coords,
	const StringT& overlap_file, double smoothing, double k2, double bound_tol, double scale_jump_0, int nip)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_22";

	/* total jump is actually 1 + scale_jump */
	scale_jump_0 -= 1.0;

	/* coarse scale element group */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	int nel = coarse->NumElements();
	if (nip == -1)
		nip = coarse->NumIP();
	else if (nip != 1)
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

	/* collect only bonds terminating with ghost points */
	RaggedArray2DT<int> ghost_neighbors_all;
	iArrayT overlap_cell_all;
	GhostNodeBonds(point_neighbors, ghost_neighbors_all, overlap_cell_all);

	/* overlap region information */
	dArray2DT bond_densities;

	/* look for restart files */
	bool solve_density = true;
	ifstreamT overlap(overlap_file);
	if (overlap.is_open())
	{
		cout << "\n " << caller << ": reading overlap from file \"" << overlap_file << '\"' << endl;
	
		/* dimension data */
		int num_cells = -99, num_bond = -99, num_ip = -99;
		overlap >> num_cells >> num_bond >> num_ip;

		/* read overlap cells indicies */
		iArrayT overlap_cell_all_tmp(num_cells);
		overlap >> overlap_cell_all_tmp;
		
		/* checks */
		if (num_bond == bonds.MajorDim() &&
		    num_ip == nip &&
		    overlap_cell_all_tmp == overlap_cell_all)
		{
			/* read bond densities */
			bond_densities.Dimension(num_cells, num_bond*num_ip);
			overlap >> bond_densities;
		
			/* don't need to solve */
			solve_density = false;
		}
		else
			cout << "\n " << caller << ": parameter mismatch in file \"" << overlap_file << '\"' << endl;

		/* close stream */
		overlap.close();
	}

	if (solve_density) /* solve for bond densities */
	{
#if 0
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
#endif

	/* set work space */
	InverseMapT overlap_cell_all_map;
	overlap_cell_all_map.SetOutOfRange(InverseMapT::MinusOne);	
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

	/* unknown bond densities */
	bond_densities.Dimension(overlap_cell_all.Length(), nip*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space that changes for each bond family */
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
	dArrayT sum_R_N;
	dArrayT f_a;
	nArrayGroupT<double> overlap_node_group(0, false);
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
	dArrayT K_scale;
	VariArrayT<double> K_scale_man(0, K_scale);
	dArrayT K_scale_last;
	VariArrayT<double> K_scale_last_man(0, K_scale);
	dArrayT diagonal;
	VariArrayT<double> diagonal_man(0, diagonal);
	dArrayT last_solution;
	VariArrayT<double> last_solution_man(0, last_solution);

	/* secant method solver */
	SecantMethodT secant_search(15, 0.05);

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
	nVariArray2DT<int> bond_densities_i_eq_man(0, bond_densities_i_eq, nip);
	RaggedArray2DT<int> inv_connects_i; 
	RaggedArray2DT<int> inv_equations_i; 
	AutoArrayT<int> bondfree_cell_i;
	bool first_pass = true;
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond family at a time */ {

		bool first_pass = true;

		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
	
		/* collect "acive" ghost node bonds - bonds of type R_i terminating at a ghost node */
		RaggedArray2DT<int> ghost_neighbors_i;
		GhostNodeBonds_2(R_i, point_coords, ghost_neighbors_all, ghost_neighbors_i, overlap_cell_i, overlap_node_i);
		overlap_cell_i_map.SetMap(overlap_cell_i);
		overlap_node_i_map.SetMap(overlap_node_i);
		
		/* report */
		cout << i+1 << '/' << bonds.MajorDim() << ": " << ghost_neighbors_i.Length() << " bonds\n";
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
		int num_eq = 0;		
		bond_densities_i_eq_man.SetMajorDimension(overlap_cell_i.Length(), false);
		bond_densities_i_eq = -1;
		for (int j = 0; j < overlap_cell_i.Length(); j++)
			for (int k = 0; k < bond_densities_i_eq.MinorDim(); k++)
				bond_densities_i_eq(j,k) = ++num_eq;

		/* "inverse" connectivities for overlap cells in the support of overlap nodes */
		TransposeConnects(*coarse, overlap_node_i, overlap_cell_i, inv_connects_i);
		
		/* collect equation numbers for "inverse" elements */
		inv_equations_i.Configure(inv_connects_i, nip);
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
		ddf_dpdp_i.AddEquationSet(inv_equations_i);
		ddf_dpdp_i.Initialize(num_eq, num_eq, 1);
		rhs_man.SetLength(num_eq, false);
		update_man.SetLength(num_eq, false);
		K_scale_man.SetLength(num_eq, false);
		K_scale_last_man.SetLength(num_eq, false);
		diagonal_man.SetLength(num_eq, false);
		last_solution_man.SetLength(num_eq, false);
		
		/* compute contribution from bonds terminating at "ghost" atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors_i, point_coords, overlap_node_i_map, sum_R_N);
		if (fPrintInput) {
			fMainOut << "ghost bond contritbution:\n";
			fMainOut << "R.sum_R_N =\n" << sum_R_N << endl;
		}

		/* initialize */
		p_i = 1.0;
		double k_penalty = 1.0e-04;
		K_scale = k_penalty;

		/* solve bond densities */
		double scale_jump = (scale_jump_0 < 1.0) ? 1.0 : scale_jump_0; /* total jump is actually 1 + scale_jump */
		double constraint_tol = bound_tol;
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+32;		
		double error_0, error;
		int max_iter = 25; /* maximum number of iterations */
		int fast_iter = 3; /* 'fast' solution to trigger bigger scale jump */
		int check_iter = 8; /* number of iterations at which the residual should have dropped to < 1.0 */
		bool more_continuation = true;
		bool do_line_search = false;
		while (more_continuation) /* continuation */		
		{		
			int iter = 0;

			/* compute residual - add Cauchy-Born contribution */
			f_a = sum_R_N;
			Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
				bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
				p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);

			/* assemble residual */
			rhs = 0.0;
			for (int j = 0; j < df_dp_i.Length(); j++)
				rhs[j] -= df_dp_i[j];

			/* copy diagonal values */
			if (!ddf_dpdp_i.CopyDiagonal(diagonal))
				ExceptionT::GeneralFail(caller);

			/* apply constraints */
			for (int j = 0; j < num_eq; j++)
			{
				double g1 = p_i[j] - 0.0;
				double g2 = 1.0 - p_i[j];
				if (g1 < 0.0)
				{
					double k = diagonal[j]*K_scale[j];
					rhs[j] -= k*g1;
					diagonal[j] = k;
				}
				else if (g2 < 0.0)
				{
					double k = diagonal[j]*K_scale[j];
					rhs[j] += k*g2;
					diagonal[j] = k;
				}
				else /* no constraint */
					diagonal[j] = 0.0;				
			}
			
			/* assemble constraints */
			ddf_dpdp_i.Assemble(diagonal, bond_densities_i_eq);

			if (false && fPrintInput) {
				fMainOut << "residual =\n" << df_dp_i << endl;
			}

			error_0 = error = sqrt(dArrayT::Dot(rhs, rhs));
			cout << i+1 << '/' << bonds.MajorDim() << ": {e, ||fa||} = {" 
			     << error << ", "
			     <<  sqrt(dArrayT::Dot(f_a, f_a)) << "}" << endl;

			while (iter < max_iter && error > abs_tol && error/error_0 > rel_tol && error/error_0 < div_tol) /* convergence */ {

				iter++;

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
					int num_active = 0;
					while (ls_continue) /* line search */
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
						s_last = s;					

						/* update densities */
						for (int j = 0; j < num_eq; j++)
							p_i[j] += rhs[j];

						/* recompute residual */			
						f_a = sum_R_N;
						Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map,
							bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
							p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);

						/* copy diagonal values */
						ddf_dpdp_i.CopyDiagonal(diagonal);
	
						/* assemble residual */
						rhs = 0.0;
						for (int j = 0; j < df_dp_i.Length(); j++)
							rhs[j] -= df_dp_i[j];
				
						/* apply constraints */
						num_active = 0;
						for (int j = 0; j < num_eq; j++)
						{
							double g1 = p_i[j] - 0.0;
							double g2 = 1.0 - p_i[j];
							if (g1 < 0.0)
							{
								num_active++;
								double k = diagonal[j]*K_scale[j];
								rhs[j] -= k*g1;
								diagonal[j] = k;
							}
							else if (g2 < 0.0)
							{
								num_active++;
								double k = diagonal[j]*K_scale[j];
								rhs[j] += k*g2;
								diagonal[j] = k;
							}
							else /* no constraint */
								diagonal[j] = 0.0;				
						}

						/* assemble constraints */
						ddf_dpdp_i.Assemble(diagonal, bond_densities_i_eq);

						/* orthogonality */
						double G_R = dArrayT::Dot(rhs, update);

						/* process results */
						if (ls_count == 0)
						{
							G_R_0 = G_R;
						
							/* no line search */
							if (!do_line_search)
								ls_continue = false;
					
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

					} /* line search */				

					error = sqrt(dArrayT::Dot(rhs, rhs));
					cout << setw(5) << iter << ": {n_c, n_ls, s, e/e_0, ||fa||} = {"
					     << setw(5) << num_active << ", " 
					     << setw(3) << secant_search.Iterations()+1 << ", " 
					     << s << ", " << error/error_0 << ", " 
					     <<  sqrt(dArrayT::Dot(f_a, f_a)) << "}" << endl;

				} /* end try */

				catch (ExceptionT::CodeT error) {
					iter = max_iter; /* exit */
				}
				
			} /* convergence */

			if (error < abs_tol || error/error_0 < rel_tol) /* converged */
			{
				/* store solution */
				K_scale_last = K_scale;
				last_solution.CopyIn(0, p_i);
			
				more_continuation = false;
				if (p_i.Length() > 0) {
					
					/* increase scaling jump */
					if (iter < fast_iter && iter > 1 && scale_jump < 200.0)
						scale_jump *= 2.0;
				
					double min = p_i[0];
					double max = p_i[0];
					for (int j = 0; j < p_i.Length(); j++) {

						/* find limits */
						double& p = p_i[j];
						min = (p < min) ? p : min;
						max = (p > max) ? p : max;
					
						/* tighter constraints */
						if (p < -constraint_tol || p - 1.0 > constraint_tol) {
							more_continuation = true;
							K_scale[j] *= (1.0 + scale_jump);
						}
					}
					
					/* report */
					cout << "{max, min, scale_jump} = " << "{" << min << ", " << max << ", " << 1+scale_jump <<"}"<< endl;
				}
			}
			else /* did not converge */
			{
				/* exit */
				if (scale_jump < 1.0e-10)
					ExceptionT::GeneralFail(caller, "scale_jump is small 1.0 + %g", scale_jump);
				else
				{
					/* restore solution */
					K_scale = K_scale_last;
					p_i.CopyIn(0, last_solution);
				
					/* reduce scaling jump */
					scale_jump *= 0.5;
				}
			}

			//cout << "{bounds, num_continuation} = {" << ", " << num_continuation << "}" << endl;

		} /* continuation */

		if (error > abs_tol && error/error_0 > rel_tol) /* failed to converged */ {
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

	/* write overlap information to restart file */
	ofstreamT overlap(overlap_file);
	overlap.precision(DBL_DIG - 1); /* full precision */
	overlap << overlap_cell_all.Length() << " "
	        << bonds.MajorDim() << " "
	        << nip << '\n';
	overlap << overlap_cell_all.wrap_tight(10) << '\n';
	overlap << bond_densities << '\n';
	overlap.close();

	} /* solve for bond densities */
	
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
		file.Root(fInputFile);
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

	/* collect current element status */
	coarse->GetStatus(fElementStatus);

	/* finding free vs projected nodes */
	int nnd = fNodeManager->NumNodes();
	ArrayT<char> node_type(nnd);
	node_type = free_;
	for (int i = 0; i < fProjectedNodes.Length(); i++) /* mark projected nodes */
		node_type[fProjectedNodes[i]] = not_free_;

	/* disable elements without any free nodes - preserve previously disabled element */
	for (int i = 0; i < nel; i++)
	{
		/* element information */
		const ElementCardT& element = coarse->ElementCard(i);
	
		/* look for free node */
		bool has_free = false;
		const iArrayT& nodes = element.NodesU();
		for (int j = 0; !has_free && j < nodes.Length(); j++)
			has_free = (node_type[nodes[j]] == free_);
		
		/* disable element */
		if (!has_free) fElementStatus[i] = ElementCardT::kOFF;
	}

	/* re-enable any overlap cells - overrides previously disabled elements */
	for (int i = 0; i < overlap_cell_all.Length(); i++)
		fElementStatus[overlap_cell_all[i]] = ElementCardT::kON;

	/* set element status */
	non_const_coarse->SetStatus(fElementStatus);
}

#endif  /* BRIDGING_ELEMENT */
