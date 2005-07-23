/* $Id: AztecParamsT.cpp,v 1.1 2005-04-05 16:07:07 paklein Exp $ */
/* created: paklein (08/12/1998) */
#include "AztecParamsT.h"

/* library support options */
#ifdef __AZTEC__

#include "az_aztec_defs.h"

using namespace Tahoe;

/* option constants and indices */
const int kNumOptions = 27; // AZ_OPTIONS_SIZE is 47, but only support some
const char* opt_name[kNumOptions] = {
"solver",
"scaling",
"preconditioner",
"convergence",
"output",
"pre_calc",
"max_iter",
"polynomial_order",
"overlap",
"type_overlap",
"k-space",
"orthogonality",
"aux_vec",
"reorder",
"keep_info",
"recursion_level",
"print_freq",
"graph_fill",
"subdomain_solve",
"init_guess",
"keep_kvecs",
"apply_kvecs",
"orth_kvecs",
"ignore_scaling",
"check_update_size",
"extreme",
"domain_solve"
};

const int opt_idex[kNumOptions] = {
AZ_solver,
AZ_scaling,
AZ_precond,
AZ_conv,
AZ_output,
AZ_pre_calc,
AZ_max_iter,
AZ_poly_ord,
AZ_overlap,
AZ_type_overlap,
AZ_kspace,
AZ_orthog,
AZ_aux_vec,
AZ_reorder,
AZ_keep_info,
AZ_recursion_level,
AZ_print_freq,
AZ_graph_fill,
AZ_subdomain_solve,
AZ_init_guess,
AZ_keep_kvecs,
AZ_apply_kvecs,
AZ_orth_kvecs,
AZ_ignore_scaling,
AZ_check_update_size,
AZ_extreme,
AZ_subdomain_solve
};

/* parameter constants and indices */
const int   kNumParams = 3; // AZ_PARAMS_SIZE is 30, but only support some
const char* param_name[kNumParams] = {"tolerance", "drop", "ilut_fill"};
const int   param_idex[kNumParams] = {AZ_tol, AZ_drop, AZ_ilut_fill};

/* constructor */
AztecParamsT::AztecParamsT(void): ParameterInterfaceT("Aztec_matrix") {}

/* set solver options */
void AztecParamsT::SetAztecOptions(int* options, double* params) const
{
	/* non-default options */
	for (int i = 0; i < AZ_options.Length(); i++)
		options[AZ_options_dex[i]] = AZ_options[i];

	/* non-default parameters */
	for (int j = 0; j < AZ_params.Length(); j++)
		params[AZ_params_dex[j]] = AZ_params[j];
}

/* describe the parameters needed by the interface */
void AztecParamsT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* solver */
	ParameterT solver(ParameterT::Enumeration, "solver");
	solver.AddEnumeration("cg", AZ_cg);
	solver.AddEnumeration("gmres", AZ_gmres);
	solver.AddEnumeration("cgs", AZ_cgs);
	solver.AddEnumeration("tfqmr", AZ_tfqmr);
	solver.AddEnumeration("bicgstab", AZ_bicgstab);
	solver.AddEnumeration("slu", AZ_slu);
	solver.AddEnumeration("symmlq", AZ_symmlq);
	solver.AddEnumeration("GMRESR", AZ_GMRESR);
	solver.AddEnumeration("fixed_pt", AZ_fixed_pt);
	solver.AddEnumeration("analyze", AZ_analyze);
	solver.AddEnumeration("lu", AZ_lu);
	solver.SetDefault(AZ_cgs);
	list.AddParameter(solver);
	
	/* scaling */
	ParameterT scaling(ParameterT::Enumeration, "scaling");
	scaling.AddEnumeration("none", AZ_none);
	scaling.AddEnumeration("Jacobi", AZ_Jacobi);
	scaling.AddEnumeration("BJacobi", AZ_BJacobi);
	scaling.AddEnumeration("row_sum", AZ_row_sum);
	scaling.AddEnumeration("sym_diag", AZ_sym_diag);
	scaling.AddEnumeration("sym_row_sum", AZ_sym_row_sum);
	scaling.AddEnumeration("equil", AZ_equil);
	scaling.AddEnumeration("sym_BJacobi", AZ_sym_BJacobi);
	scaling.SetDefault(AZ_none);
	list.AddParameter(scaling);

	/* preconditioner */
	ParameterT precond(ParameterT::Enumeration, "preconditioner");
	precond.AddEnumeration("none", AZ_none);
	precond.AddEnumeration("Jacobi", AZ_Jacobi);
	precond.AddEnumeration("sym_GS", AZ_sym_GS);
	precond.AddEnumeration("Neumann", AZ_Neumann);
	precond.AddEnumeration("ls", AZ_ls);
	precond.AddEnumeration("ilu", AZ_ilu);
	precond.AddEnumeration("bilu", AZ_bilu);
	precond.AddEnumeration("icc", AZ_icc);
	precond.AddEnumeration("ilut", AZ_ilut);
	precond.AddEnumeration("rilu", AZ_rilu);
	precond.AddEnumeration("recursive", AZ_recursive);
	precond.AddEnumeration("smoother", AZ_smoother);
	precond.AddEnumeration("dom_decomp", AZ_dom_decomp);
	precond.AddEnumeration("user_precond", AZ_user_precond);
	precond.AddEnumeration("bilu_ifp", AZ_bilu_ifp);
	precond.SetDefault(AZ_ls);
	list.AddParameter(precond);

	/* convergence */
	ParameterT conv(ParameterT::Enumeration, "convergence");
	conv.AddEnumeration("r0", AZ_r0);
	conv.AddEnumeration("rhs", AZ_rhs);
	conv.AddEnumeration("Anorm", AZ_Anorm);
	conv.AddEnumeration("sol", AZ_sol);
	conv.AddEnumeration("weighted", AZ_weighted);
	conv.AddEnumeration("expected_values", AZ_expected_values);
	conv.AddEnumeration("noscaled", AZ_noscaled);
	conv.AddEnumeration("AZTECOO_conv_test", AZTECOO_conv_test);
	conv.AddEnumeration("inf_noscaled", AZ_inf_noscaled);
	conv.SetDefault(AZ_r0);
	list.AddParameter(conv);

	/* output */
	ParameterT output(ParameterT::Enumeration, "output");
	output.AddEnumeration("1", 1);
	output.AddEnumeration("5", 5);
	output.AddEnumeration("10", 10);
	output.AddEnumeration("25", 25);
	output.AddEnumeration("50", 50);
	output.AddEnumeration("100", 100);
	output.AddEnumeration("500", 500);
	output.AddEnumeration("1000", 1000);
	output.AddEnumeration("none", AZ_none);
	output.AddEnumeration("all", AZ_all);
	output.AddEnumeration("last", AZ_last);
	output.AddEnumeration("warnings", AZ_warnings);
	output.SetDefault(AZ_none);
	list.AddParameter(output);

	/* factorization */
	ParameterT pre_calc(ParameterT::Enumeration, "pre_calc");
	pre_calc.AddEnumeration("calc", AZ_calc);
	pre_calc.AddEnumeration("recalc", AZ_recalc);
	pre_calc.AddEnumeration("reuse", AZ_reuse);
	pre_calc.AddEnumeration("sys_reuse", AZ_sys_reuse);
	pre_calc.SetDefault(AZ_calc);
	list.AddParameter(pre_calc);
	
	ParameterT max_iter(ParameterT::Integer, "max_iter");
	max_iter.AddLimit(0, LimitT::Lower);
	max_iter.SetDefault(500);
	list.AddParameter(max_iter);
	
	ParameterT poly_ord(ParameterT::Integer, "polynomial_order");
	poly_ord.SetDefault(5);
	list.AddParameter(poly_ord);

	/* domain decomposition overlap */
	ParameterT overlap(ParameterT::Enumeration, "overlap");
	overlap.AddEnumeration("none", AZ_none);
	overlap.AddEnumeration("diag", AZ_diag);
	overlap.AddEnumeration("full", AZ_full);
	overlap.SetDefault(AZ_none);
	list.AddParameter(overlap);

	ParameterT kspace(ParameterT::Integer, "k-space");
	kspace.SetDefault(30);
	list.AddParameter(kspace);

	ParameterT orthog(ParameterT::Enumeration, "orthogonality");
	orthog.AddEnumeration("classic", AZ_classic);
	orthog.AddEnumeration("modified", AZ_modified);
	orthog.SetDefault(AZ_modified);
	list.AddParameter(orthog);

	ParameterT aux_vec(ParameterT::Enumeration, "aux_vec");
	aux_vec.AddEnumeration("resid", AZ_resid);
	aux_vec.AddEnumeration("rand", AZ_rand);
	aux_vec.SetDefault(AZ_resid);
	list.AddParameter(aux_vec);

	ParameterT tol(ParameterT::Double, "tolerance");
	tol.SetDefault(4.00e-9);
	list.AddParameter(tol);

	ParameterT drop(ParameterT::Double, "drop");
	drop.SetDefault(0.0);
	list.AddParameter(drop);

	ParameterT ilut_fill(ParameterT::Double, "ilut_fill");
	ilut_fill.SetDefault(0.0);
	list.AddParameter(ilut_fill);

	ParameterT subdomain_solve(ParameterT::Integer, "subdomain_solve");
	list.AddParameter(subdomain_solve);
}

/* describe the parameters needed by the interface */
void AztecParamsT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	AutoArrayT<int> opt, opt_dex, param_dex;
	AutoArrayT<double> param;
	const ArrayT<ParameterT>& parameters = list.Parameters();
	for (int i = 0; i < parameters.Length(); i++) {
		const ParameterT& parameter = parameters[i];
		int o_dex = OptionNameToIndex(parameter.Name());
		int p_dex = ParamNameToIndex(parameter.Name());
		if (o_dex != -1) {
			int o = parameter;
			opt.Append(o);
			opt_dex.Append(o_dex);
		}
		else if (p_dex != -1) {
			double p = parameter;
			param.Append(p);
			param_dex.Append(p_dex);
		}
		else
			ExceptionT::GeneralFail("AztecParamsT::TakeParameterList",
				"could not resolve \"%s\"", parameter.Name().Pointer());
	}
	
	/* keep */
	AZ_options_dex.Swap(opt_dex);
	AZ_options.Swap(opt);
	AZ_params_dex.Swap(param_dex);
	AZ_params.Swap(param);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* option name to index conversion */
int AztecParamsT::OptionNameToIndex(const char* name) const
{
	int i, found = 0;
	for (i = 0; i < kNumOptions && !found; i++)
		if (strcmp(name,opt_name[i]) == 0) found = 1;

	if (found)
		return opt_idex[i-1];
	else
		return -1;
}

/* parameter name to index conversion */
int AztecParamsT::ParamNameToIndex(const char* name) const
{
	int i, found = 0;
	for (i = 0; i < kNumParams && !found; i++)
		if (strcmp(name,param_name[i]) == 0) found = 1;

	if (found)
		return param_idex[i-1];
	else
		return -1;
}

/* library support options */
#endif /* __AZTEC__ */
