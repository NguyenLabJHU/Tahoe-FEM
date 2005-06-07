/* $Id: FEManagerT_bridging_2.cpp,v 1.14 2005-04-28 23:57:15 paklein Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
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

using namespace Tahoe;

void FEManagerT_bridging::CorrectOverlap_2(const RaggedArray2DT<int>& point_neighbors, const dArray2DT& point_coords,
	const StringT& overlap_file, double smoothing, double k2, double k_r, double bound_0, int nip)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_2";

	//TEMP - currently unused
	int do_line_search = 0;

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
	if (fLogging == GlobalT::kVerbose) {
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
	SPOOLESMatrixT ddf_dpdp_i(Output(), GlobalMatrixT::kZeroPivots, true, true, 0, fComm);
#else
	CCSMatrixT ddf_dpdp_i(Output(), GlobalMatrixT::kZeroPivots, fComm);
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
	dArrayT last_solution;
	VariArrayT<double> last_solution_man(0, last_solution);
	dArrayT constraint;
	VariArrayT<double> constraint_man(0, constraint);

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
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond family at a time */ {

		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
	
		/* collect "acive" ghost node bonds - bonds of type R_i terminating at a ghost node */
		RaggedArray2DT<int> ghost_neighbors_i;
		GhostNodeBonds_2(R_i, point_coords, ghost_neighbors_all, ghost_neighbors_i, overlap_cell_i, overlap_node_i);
		overlap_cell_i_map.SetMap(overlap_cell_i);
		overlap_node_i_map.SetMap(overlap_node_i);

		/* report */
		cout << i+1 << '/' << bonds.MajorDim() << ": " << ghost_neighbors_i.Length() << " bonds\n";
		if (fLogging == GlobalT::kVerbose) {
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
		else if (fLogging != GlobalT::kSilent)
			fMainOut << i+1 << '/' << bonds.MajorDim() << ": " << ghost_neighbors_i.Length() << " bonds\n";
		
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
		bond_densities_i_eq_man.SetMajorDimension(overlap_cell_i.Length(), false);
		bond_densities_i_eq = -1;
		int num_eq = 0;
		for (int j = 0; j < overlap_cell_i.Length(); j++)
			for (int k = 0; k < bond_densities_i_eq.MinorDim(); k++)
				bond_densities_i_eq(j,k) = ++num_eq;
		int num_eq_p = num_eq;				

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

		/* configure linear solver */
		ddf_dpdp_i.AddEquationSet(inv_equations_i);
		ddf_dpdp_i.Initialize(num_eq, num_eq, 1);
		rhs_man.SetLength(num_eq, false);
		update_man.SetLength(num_eq, false);
		last_solution_man.SetLength(num_eq, false);
		
		/* compute contribution from bonds terminating at "ghost" atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors_i, point_coords, overlap_node_i_map, sum_R_N);
		if (fLogging == GlobalT::kVerbose) {
			fMainOut << "ghost bond contribution:\n";
			fMainOut << "R.sum_R_N =\n" << sum_R_N << endl;
		}

		/* initialize */
		p_i = 1.0;

		/* initial number of continuation steps */
		double bound = (k_r > kSmall) ? bound_0 : 0.5;
		bool set_outer_bounds = false;
		double outer_bound = bound;
		double last_bound = bound;
		int num_continuation_0 = 2;
		int num_continuation = num_continuation_0;

		/* solve bond densities */
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+03;		
		double error_0, error;
		int max_iter = 15;
		int check_iter = 8;
		bool more_continuation = true;
		while (more_continuation) /* continuation */		
		//while (bound - 0.50 > kSmall) /* continuation */
		{		
			int iter = 0;

			/* compute residual - add Cauchy-Born contribution */
			f_a = sum_R_N;
			Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
				bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
				p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);

//DEBUGGING - check df_dp_i using finite difference
#if 0
int prec = fMainOut.precision();
fMainOut.precision(12);
fMainOut << "(1) df_dp =\n" << df_dp_i << '\n';
double p_eps = 1.0e-08;
double fa2 = dArrayT::Dot(f_a, f_a);
dArray2DT df_dp_i_tmp = df_dp_i;
for (int ii = 0; ii < p_i.Length(); ii++)
{
	/* compute overlap */
	p_i[ii] += p_eps;
	f_a = sum_R_N;
	Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
		bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
		p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
	p_i[ii] -= p_eps;

	/* compute */
	df_dp_i_tmp[ii] = (dArrayT::Dot(f_a, f_a) - fa2)/p_eps;
}
fMainOut << "(2) df_dp =\n" << df_dp_i << '\n';

/* restore residual */
f_a = sum_R_N;
Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map, 
	bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
	p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
fMainOut.precision(prec);
#endif
//DEBUGGING - check df_dp_i using finite difference

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
				double   g = bound*bound - (p_i[j] - 0.5)*(p_i[j] - 0.5);
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

			if (fLogging == GlobalT::kVerbose) {
				fMainOut << "residual =\n" << df_dp_i << endl;
			}

			error_0 = error = sqrt(dArrayT::Dot(rhs, rhs));
			cout << i+1 << '/' << bonds.MajorDim() << ": {e, bound, ||fa||} = {" 
			     << error << ", " 
			     << bound << ", " 
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
						for (int j = 0; j < num_eq_p; j++)
							p_i[j] += rhs[j];
								
						/* update constraints */
						for (int j = 0; j < num_eq_L; j++)
							constraint[j] += rhs[j+num_eq_p];

						/* recompute residual */			
						f_a = sum_R_N;
						Compute_df_dp_2(R_i, V_0, cell_type, overlap_cell_i_map, overlap_node_i, overlap_node_i_map,
							bond_densities_i_eq, inv_connects_i, inv_equations_i, inv_equations_i,
							p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);

						/* assemble residual */
						rhs = 0.0;
						for (int j = 0; j < df_dp_i.Length(); j++)
							rhs[j] -= df_dp_i[j];
				
						/* apply constraints */
						num_active = 0;
						for (int j = 0; j < num_eq_L; j++)
						{
							/* equation numbers */
							constraint_eq.RowAlias(j, eqnos);

							/* constraint equation: g >= 0 */
							double   g = bound*bound - (p_i[j] - 0.5)*(p_i[j] - 0.5);
							double  Dg = -2.0*(p_i[j] - 0.5);
							double DDg = -2.0;
		
							/* augmented multiplier */
							double L_r = constraint[j] + k_r*g;
							if (L_r < 0.0) /* active constraint */
							{
								/* count active constraints */
								num_active++;
						
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
				
				/* check progress */
				if (iter > check_iter && error/error_0 > 1.0)
					iter = max_iter; /* exit */
				
			} /* convergence */

			if (error < abs_tol || error/error_0 < rel_tol) /* converged */
			{
				/* converged at least once */
				set_outer_bounds = true;

				/* store the solution */
				last_bound = bound;
				last_solution.CopyIn(0, p_i);
				last_solution.CopyIn(num_eq_p, constraint);

				/* increase step size */
				if (iter > 1 && iter < 4 && num_continuation > num_continuation_0)
					num_continuation /= 2;

				/* tighten bounds */
				if (bound > 0.5)
				{
					double d_bound = (outer_bound - 0.5)/num_continuation;
					bound -= d_bound;
					
					
					/* correct overshoot */
					bound = (bound < 0.5) ? 0.5 : bound;
				}
				else /* done */
					more_continuation = false;
			}
			else /* did not converge */
			{
				if (!set_outer_bounds) /* expanding bounds */
				{
					bound *= 2.0;
					outer_bound = bound;
					
					/* reset solution */
					p_i = 1.0;
					constraint = 0.0;
				}
				else /* add more continuation steps */
				{
					/* reset the solution */
					bound = last_bound;
					p_i.CopyPart(0, last_solution, 0, num_eq_p);				
					constraint.CopyPart(0, last_solution, num_eq_p, num_eq_L);
				
					/* more continuation steps */
					num_continuation *= 2;
				}
			}

			cout << "{bounds, num_continuation} = {" << bound  << ", " << num_continuation << "}" << endl;

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
	if (fLogging != GlobalT::kSilent) {
	
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

void FEManagerT_bridging::GhostNodeBonds_2(const dArrayT& R_i, const dArray2DT& point_coords, 
	const RaggedArray2DT<int>& ghost_neighbors_all, RaggedArray2DT<int>& ghost_neighbors_i, 
	AutoArrayT<int>& overlap_cell_i, AutoArrayT<int>& overlap_node_i) const
{
	const char caller[] = "FEManagerT_bridging::GhostNodeBonds_2";

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

//NOTE: up to this point, overlap cells are classified as those containing an active ghost
//      atom bond; while overlap nodes are those nodes whose support includes a ghost atom.
//      Now, we grow the list of overlap cells to include the entire support of all nodes
//      whose support contains an active ghost node bond. This must be consistent with the
//      other version of FEManagerT_bridging::GhostNodeBonds used to collect the entire,
//      potential overlap region.
#if 0
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
	overlap_cell_i.Dimension(num_overlap);
	num_overlap = 0;
	for (int i = 0; i < is_overlap_cell.Length(); i++)
		if (is_overlap_cell[i] == 't')
			overlap_cell_i[num_overlap++] = i;
#endif
}

void FEManagerT_bridging::Compute_df_dp_2(const dArrayT& R, double V_0, const ArrayT<char>& cell_type, 
	const InverseMapT& overlap_cell_map, const ArrayT<int>& overlap_node, const InverseMapT& overlap_node_map,
	const iArray2DT& cell_eq_active_i,
	const RaggedArray2DT<int>& inv_connects_i, const RaggedArray2DT<int>& inv_equations_all_i, const RaggedArray2DT<int>& inv_equations_active_i,
	const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, GlobalMatrixT& ddf_dpdp) const
{
	const char caller[] = "FEManagerT_bridging::Compute_df_dp_2";

	/* the continuum element solving the coarse scale */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	if (!coarse) ExceptionT::GeneralFail(caller, "interpolation data not set");

	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse->NumElementNodes();
	int nip = rho.MinorDim();

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse->GeometryCode(), nip, element_coords);
	shapes.Initialize();

	/* B_hat_U_U */
	InterpolationDataT& B_hatU_U = const_cast<InterpolationDataT&>(fDrivenCellData.NodeToNode());
	const RaggedArray2DT<int>& B_hatU_U_neighbors = B_hatU_U.Neighbors();	
	const RaggedArray2DT<double>& B_hatU_U_weights = B_hatU_U.NeighborWeights();	
	InverseMapT& B_hatU_U_row_map = B_hatU_U.Map();
	InverseMapT::SettingT old_out_of_range = B_hatU_U_row_map.OutOfRange();
	B_hatU_U_row_map.SetOutOfRange(InverseMapT::MinusOne); /* need this to differentiate free/prescribed nodes */
	iArrayT hatU_U_neighbors;
	dArrayT hatU_U_weights;

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

	/* add B_hatU_U to free nodes from prescribed nodes */
	if (B_hatU_U_weights.Length() > 0) /* have hatU-U contributions */
	{
		for (int i = 0; i < overlap_node.Length(); i++)
		{
			int hatU_node = overlap_node[i];
			int hatU_row = B_hatU_U_row_map.Map(hatU_node);
			if (hatU_row != -1) /* is prescribed */
			{
				/* row of B_hatU_U */
				B_hatU_U_neighbors.RowAlias(hatU_row, hatU_U_neighbors);
				B_hatU_U_weights.RowAlias(hatU_row, hatU_U_weights);
				
				/* add contributions to U nodes */
				for (int j = 0; j < hatU_U_neighbors.Length(); j++)
				{
					int U_node = hatU_U_neighbors[j];
					int U_node_index = overlap_node_map.Map(U_node);
					if (U_node_index > -1) /* allowed to be -1? */
						f_a[U_node_index] += f_a[i]*hatU_U_weights[j]; /* add B_hatU_U contribution */
				}
			}
		}
	}

	/* output */
	if (0 && fLogging == GlobalT::kVerbose) fMainOut << "f_a =\n" << f_a << '\n';

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
			cell_eq_active_i.RowAlias(overlap_cell_index, eqnos);
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
					local_node = k;
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

		/* output */
		if (0 && fLogging == GlobalT::kVerbose) fMainOut << "df_a_dp =\n" << df_a_dp << '\n';

		/* assemble stiffness */
		df_a_dp_2_man.SetDimensions(inv_equations_active_i.MinorDim(i));
		df_a_dp_2.Outer(df_a_dp, df_a_dp);
		inv_equations_active_i.RowAlias(i, eqnos);
		ddf_dpdp.Assemble(df_a_dp_2, eqnos);
		
		/* force contribution */
		inv_equations_all_i.RowAlias(i, eqnos);
		for (int j = 0; j < eqnos.Length(); j++)
			df_dp[eqnos[j]-1] += f_a[node_index]*df_a_dp[j];

		if (B_hatU_U_weights.Length() > 0) /* have hatU-U contributions */
		{
			int hatU_row = B_hatU_U_row_map.Map(node);
			if (hatU_row != -1) /* is prescribed */
			{
				/* row of B_hatU_U */
				B_hatU_U_neighbors.RowAlias(hatU_row, hatU_U_neighbors);
				B_hatU_U_weights.RowAlias(hatU_row, hatU_U_weights);
			
				/* sum over U nodes */
				double Bxf = 0.0;
				for (int k = 0; k < hatU_U_neighbors.Length(); k++) {
					int U_node_index = overlap_node_map.Map(hatU_U_neighbors[k]);
					if (U_node_index > -1) /* allowed to be -1? */
						Bxf += hatU_U_weights[k]*f_a[U_node_index];
				}
			
				/* assemble */
				for (int j = 0; j < eqnos.Length(); j++)
					df_dp[eqnos[j]-1] += Bxf*df_a_dp[j];
			}
		}
	}

	/* output */
	if (0 && fLogging == GlobalT::kVerbose) fMainOut << "df_dp =\n" << df_dp << '\n';

	/* restore map behavior */
	B_hatU_U_row_map.SetOutOfRange(old_out_of_range);
}

#endif  /* BRIDGING_ELEMENT */
