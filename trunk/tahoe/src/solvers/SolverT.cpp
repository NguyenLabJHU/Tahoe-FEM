/* $Id: SolverT.cpp,v 1.22 2004-07-15 08:31:51 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "SolverT.h"

#include <iostream.h>
#include <string.h>

#include "ofstreamT.h"
#include "FEManagerT.h"
#include "CommunicatorT.h"
#include "iArrayT.h"
#include "ElementMatrixT.h"
#include "ParameterContainerT.h"

/* global matrix */
#include "CCSMatrixT.h"
#include "DiagonalMatrixT.h"
#include "FullMatrixT.h"
#include "CCNSMatrixT.h"
#include "AztecMatrixT.h"
#include "SuperLUMatrixT.h"
#include "SuperLU_DISTMatrixT.h"
#include "SPOOLESMatrixT.h"
#include "PSPASESMatrixT.h"

#ifdef __TAHOE_MPI__
#include "SPOOLESMatrixT_mpi.h"
#endif

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<SolverT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<SolverT*>::fByteCopy = true;
} /* namespace Tahoe */

SolverT::SolverT(FEManagerT& fe_manager, int group):
	ParameterInterfaceT("solver"),
	fFEManager(fe_manager),
	fGroup(group),
	fLHS(NULL),
	fNumIteration(0),
	fLHS_lock(kOpen),
	fLHS_update(true),
	fRHS_lock(kOpen),
	fPerturbation(0.0)
{
	/* console */
	iSetName("solver");
	iAddVariable("print_equation_numbers", fPrintEquationNumbers);
}

/* destructor */
SolverT::~SolverT(void) { delete fLHS; }

/* configure the global equation system */
void SolverT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{	
	try {
		/* allocate rhs vector */
		fRHS.Dimension(loc_num_eq);
		fRHS = 0.0;
		
		/* set global equation matrix type */
		fLHS->Initialize(tot_num_eq, loc_num_eq, start_eq);
	
		/* output global equation number for each DOF */
		if (fPrintEquationNumbers) fFEManager.WriteEquationNumbers(fGroup);
	}	

	catch (ExceptionT::CodeT error_code)
	{
		cout << "\n SolverT::Initialize: exception: "
		     << ExceptionT::ToString(error_code) << endl;
		throw error_code;
	}
}

/* start solution step */
void SolverT::InitStep(void)
{
	fNumIteration = -1;
	fLHS_update = true;
}

/* end solution step */
void SolverT::CloseStep(void)
{
	/* do nothing */ 
}

/* error handler */
void SolverT::ResetStep(void) 
{
	/* do nothing */ 
}

/* process element group equation data to configure GlobalMatrixT */
void SolverT::ReceiveEqns(const iArray2DT& equations) const
{
	fLHS->AddEquationSet(equations);
}

void SolverT::ReceiveEqns(const RaggedArray2DT<int>& equations) const
{
	fLHS->AddEquationSet(equations);
}

void SolverT::AssembleRHS(const nArrayT<double>& elRes, const nArrayT<int>& eqnos)
{
	const char caller[] = "SolverT::AssembleRHS";

	/* consistency check */
	if (elRes.Length() != eqnos.Length()) ExceptionT::GeneralFail(caller);

	/* lock state */
	if (fRHS_lock == kIgnore)
		return;
	else if (fRHS_lock == kLocked)
		ExceptionT::GeneralFail(caller, "RHS is locked");

#if __option(extended_errorcheck)
	GlobalT::EquationNumberScopeT scope = (GlobalT::EquationNumberScopeT) fLHS->EquationNumberScope();
#endif

	int num_eq = fLHS->NumEquations();
	int start_eq = fLHS->StartEquation();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = eqnos[i] - start_eq;
	
		/* in range */
		if (eq > -1 && eq < num_eq) fRHS[eq] += elRes[i];

#if __option(extended_errorcheck)
		else if (scope == GlobalT::kLocal && eq >= num_eq)
			ExceptionT::OutOfRange(caller, "equation number is out of range: %d", eq + start_eq);
#endif
	}	
}

void SolverT::OverWriteRHS(const dArrayT& elRes, const nArrayT<int>& eqnos)
{
	/* consistency check */
	if (elRes.Length() != eqnos.Length()) throw ExceptionT::kGeneralFail;

	/* lock state */
	if (fRHS_lock == kIgnore)
		return;
	else if (fRHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::OverWriteRHS", "RHS is locked");

	int num_eq = fLHS->NumEquations();
	int start_eq = fLHS->StartEquation();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = eqnos[i] - start_eq;
	
		/* in range */
		if (eq > -1 && eq < num_eq) fRHS[eq] = elRes[i];
	}	
}

void SolverT::DisassembleRHS(dArrayT& elRes, const nArrayT<int>& eqnos) const
{
	/* consistency check */
	if (elRes.Length() != eqnos.Length()) throw ExceptionT::kGeneralFail;

	int num_eq = fLHS->NumEquations();
	int start_eq = fLHS->StartEquation();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eq = eqnos[i] - start_eq;
	
		/* in range */
		if (eq > 0 && eq < num_eq)
			elRes[i] = fRHS[eq];
		else
			elRes[i] = 0.0;
	}	
}

/* return the required equation numbering scope - local by default */
GlobalT::EquationNumberScopeT SolverT::EquationNumberScope(void) const
{
#if __option(extended_errorcheck)
	if (!fLHS)
	{
		cout << "\n SolverT::EquationNumberScope: invalid LHS" << endl;
		throw ExceptionT::kGeneralFail;	
	}
#endif

	return (GlobalT::EquationNumberScopeT) fLHS->EquationNumberScope();
}

/* describe the parameters needed by the interface */
void SolverT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* print equation numbers */
	ParameterT print_eqnos(ParameterT::Boolean, "print_eqnos");
	print_eqnos.SetDefault(false);
	list.AddParameter(print_eqnos);

	/* check code */
	ParameterT check_code(ParameterT::Enumeration, "check_code");
	check_code.AddEnumeration("no_check", GlobalMatrixT::kNoCheck);
	check_code.AddEnumeration("small_pivots", GlobalMatrixT::kZeroPivots);
	check_code.AddEnumeration("print_LHS", GlobalMatrixT::kPrintLHS);
	check_code.AddEnumeration("print_RHS", GlobalMatrixT::kPrintRHS);
	check_code.AddEnumeration("print_solution", GlobalMatrixT::kPrintSolution);
	check_code.AddEnumeration("check_LHS", GlobalMatrixT::kCheckLHS);
	check_code.SetDefault(GlobalMatrixT::kNoCheck);
	list.AddParameter(check_code);
	
	/* perturbation used to compute LHS check */
	ParameterT check_LHS_perturbation(fPerturbation, "check_LHS_perturbation");
	check_LHS_perturbation.AddLimit(0.0, LimitT::LowerInclusive);
	check_LHS_perturbation.SetDefault(1.0e-08);
	list.AddParameter(check_LHS_perturbation);
}

/* information about subordinate parameter lists */
void SolverT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* linear solver choice */
	sub_list.AddSub("matrix_type_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SolverT::NewSub(const StringT& name) const
{
	if (name == "matrix_type_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
	
		choice->AddSub(ParameterContainerT("profile_matrix"));
		choice->AddSub(ParameterContainerT("diagonal_matrix"));
		choice->AddSub(ParameterContainerT("full_matrix"));

#ifdef __AZTEC__
		choice->AddSub(ParameterContainerT("Aztec_matrix"));
#endif

#ifdef __SPOOLES__
		choice->AddSub(ParameterContainerT("SPOOLES_matrix"));
#endif

#ifdef __PSPASES__
		choice->AddSub(ParameterContainerT("PSPASES_matrix"));
#endif	

#ifdef __SUPERLU__
		choice->AddSub(ParameterContainerT("SuperLU_matrix"));
#endif	

		return choice;
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void SolverT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);
	
	/* construct matrix */
	fPrintEquationNumbers = list.GetParameter("print_eqnos");
	int check_code = list.GetParameter("check_code");
	fPerturbation = list.GetParameter("check_LHS_perturbation");
	SetGlobalMatrix(list.GetListChoice(*this, "matrix_type_choice"), check_code);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* return the magnitude of the residual force */
double SolverT::Residual(const dArrayT& force) const
{
	return sqrt(InnerProduct(force,force));
}

/* (distributed) inner product */	
double SolverT::InnerProduct(const dArrayT& v1, const dArrayT& v2) const
{
	/* check heart beat */
	if (fFEManager.Communicator().Sum(ExceptionT::kNoError) != 0) throw ExceptionT::kBadHeartBeat;

	return fFEManager.Communicator().Sum(dArrayT::Dot(v1, v2));
}

/* return approximate stiffness matrix */
GlobalMatrixT* SolverT::ApproximateLHS(const GlobalMatrixT& template_LHS)
{
	/* create matrix with same structure as the template */
	GlobalMatrixT* approx_LHS = template_LHS.Clone();

	/* open locks */
	fRHS_lock = kOpen;
	fLHS_lock = kIgnore;

	/* get copy of residual */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	dArrayT rhs = fRHS;
	dArrayT update;
	update.Dimension(rhs);
	update = 0.0;
	
	/* perturb each degree of freedom and compute the new residual */
	approx_LHS->Clear();
	iArrayT col(1);
	AutoArrayT<int> rows;
	AutoArrayT<double> K_col_tmp;
	ElementMatrixT K_col(ElementMatrixT::kNonSymmetric);
	for (int i = 0; i < fRHS.Length(); i++)
	{
		/* perturbation */
		update[i] = fPerturbation;
		
		/* apply update to system */
		fFEManager.Update(Group(), update);

		/* compute residual */
		fRHS = 0.0;
		fFEManager.FormRHS(Group());	
	
		/* reset work space */
		rows.Dimension(0);
		K_col_tmp.Dimension(0);
			
		/* compute column of stiffness matrix */
		for (int j = 0; j < fRHS.Length(); j++)
		{
			/* finite difference approximation */
			double K_ij = (rhs[j] - fRHS[j])/fPerturbation;

			/* assemble only non-zero values */
			if (fabs(K_ij) > kSmall) {
				col[0] = i+1;
				rows.Append(j+1);
				K_col_tmp.Append(K_ij);
			}
		}
		
		/* assemble */
		K_col.Alias(rows.Length(), 1, K_col_tmp.Pointer());
		approx_LHS->Assemble(K_col, rows, col);
			
		/* undo perturbation */
		update[i] = -fPerturbation;
		if (i > 0) update[i-1] = 0.0;
	}
		
	/* restore configuration and residual */
	fFEManager.Update(Group(), update);
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	

	/* close locks */
	fRHS_lock = kLocked;
	fLHS_lock = kLocked;

	/* return */
	return approx_LHS;
}

void SolverT::CompareLHS(const GlobalMatrixT& ref_LHS, const GlobalMatrixT& test_LHS) const
{
	ofstreamT& out = fFEManager.Output();

	out << "\nreference LHS:\n";
	ref_LHS.PrintLHS(true);

	out << "\ntest LHS:\n";
	test_LHS.PrintLHS(true);
}

/*************************************************************************
* Private
*************************************************************************/

/* check matrix type against analysis code, return 1 if
* compatible, 0 otherwise */
int SolverT::CheckMatrixType(int matrix_type, int analysis_code) const
{
	int OK = 1;
	switch (matrix_type)
	{
		case kDiagonalMatrix:
		
			OK = (analysis_code == GlobalT::kLinExpDynamic   ||
			      analysis_code == GlobalT::kNLExpDynamic    ||
			      analysis_code == GlobalT::kVarNodeNLExpDyn ||
			      analysis_code == GlobalT::kNLExpDynKfield);
			break;
		
		case kProfileSolver:
		
			OK = (analysis_code == GlobalT::kLinStatic       ||
			      analysis_code == GlobalT::kLinDynamic      ||
			      analysis_code == GlobalT::kNLStatic        ||
			      analysis_code == GlobalT::kNLDynamic       ||
			      analysis_code == GlobalT::kDR              ||
			      analysis_code == GlobalT::kNLStaticKfield  ||
			      analysis_code == GlobalT::kVarNodeNLStatic ||
			      analysis_code == GlobalT::kAugLagStatic);
			break;

		case kFullMatrix:

			/* not for explicit dynamics */
			OK = (analysis_code != GlobalT::kLinExpDynamic &&
			      analysis_code != GlobalT::kNLExpDynamic  &&
			      analysis_code != GlobalT::kVarNodeNLExpDyn);
			break;
					
		case kAztec:

			/* not for explicit dynamics */
			OK = (analysis_code != GlobalT::kLinExpDynamic &&
			      analysis_code != GlobalT::kNLExpDynamic  &&
			      analysis_code != GlobalT::kVarNodeNLExpDyn);
			break;
			
		case kSuperLU:
		
			OK = (analysis_code == GlobalT::kLinStatic       ||
			      analysis_code == GlobalT::kLinDynamic      ||
			      analysis_code == GlobalT::kNLStatic        ||
			      analysis_code == GlobalT::kNLDynamic       ||
			      analysis_code == GlobalT::kDR              ||
			      analysis_code == GlobalT::kNLStaticKfield  ||
			      analysis_code == GlobalT::kVarNodeNLStatic ||
			      analysis_code == GlobalT::kAugLagStatic);
			break;

		case kSPOOLES:
		
			OK = (analysis_code == GlobalT::kLinStatic       ||
			      analysis_code == GlobalT::kLinDynamic      ||
			      analysis_code == GlobalT::kNLStatic        ||
			      analysis_code == GlobalT::kNLDynamic       ||
			      analysis_code == GlobalT::kDR              ||
			      analysis_code == GlobalT::kNLStaticKfield  ||
			      analysis_code == GlobalT::kVarNodeNLStatic ||
			      analysis_code == GlobalT::kAugLagStatic);		
			break;

		default:

			cout << "\n SolverT::CheckMatrixType: unknown matrix type ";
			cout << matrix_type << '\n';
			OK = 0;		
	}

	/* compatibility */
	if (!OK)
	{
		cout << "\n SolverT::CheckMatrixType: matrix type " << matrix_type << '\n';
		cout << " is not compatible with analysis code " << analysis_code << endl;

		ostream& out = fFEManager.Output();
		out << "\n SolverT::CheckMatrixType: matrix type " << matrix_type << '\n';
		out << " is not compatible with analysis code " << analysis_code << endl;
	}
	else
		return 1;

	/* checks */
	if (fFEManager.Size() > 1 &&
	    (matrix_type == kFullMatrix    ||
	     matrix_type == kProfileSolver ||
	     matrix_type == kSuperLU  ||
	     matrix_type == kSPOOLES))
	{
		cout << "\n SolverT::CheckMatrixType: matrix type not support in parallel: "
		     << matrix_type << endl;
		throw ExceptionT::kGeneralFail;
	}
	return OK;
}

/* set global equation matrix */
void SolverT::SetGlobalMatrix(const ParameterListT& params, int check_code)
{
	const char caller[] = "SolverT::SetGlobalMatrix";

	/* streams */
//	ifstreamT&  in = fFEManager.Input();
	ofstreamT& out = fFEManager.Output();

	/* resolve matrix type */
	if (params.Name() == "diagonal_matrix")
	{
		DiagonalMatrixT* diag = new DiagonalMatrixT(out, check_code, DiagonalMatrixT::kNoAssembly);
		diag->SetAssemblyMode(DiagonalMatrixT::kDiagOnly);
		fLHS = diag;
	}
	else if (params.Name() == "profile_matrix")
	{
		/* global system properties */
		GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);
		
		if (type == GlobalT::kNonSymmetric)
			fLHS = new CCNSMatrixT(out, check_code);
		else if (type == GlobalT::kSymmetric)
			fLHS = new CCSMatrixT(out, check_code);
		else
			ExceptionT::GeneralFail(caller, "system type %d not compatible with matrix %d", type, kProfileSolver);
	}
	else if (params.Name() == "full_matrix")
		fLHS = new FullMatrixT(out, check_code);
	else if (params.Name() == "SPOOLES_matrix")
	{
#ifdef __SPOOLES__
		/* global system properties */
		GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);

		/* solver options */
		bool pivoting = true; //NOTE: SPOOLES v2.2 does not seem to solve non-symmetric
			                      //      systems correctly in parallel if pivoting is disabled
		bool symmetric;
		if (type == GlobalT::kDiagonal || type == GlobalT::kSymmetric)
			symmetric = true;
		else if (type == GlobalT::kNonSymmetric)
			symmetric = false;
		else
			ExceptionT::GeneralFail(caller, "unexpected system type: %d", type);

#ifdef __TAHOE_MPI__
		/* constuctor */
		if (fFEManager.Size() > 1)
		{
#ifdef __SPOOLES_MPI__
			fLHS = new SPOOLESMatrixT_mpi(out, check_code, symmetric, pivoting, fFEManager.Communicator());
#else /* __SPOOLES_MPI__ */
			ExceptionT::GeneralFail(caller, "SPOOLES MPI not installed");
#endif /* __SPOOLES_MPI__ */
		}
		else
			fLHS = new SPOOLESMatrixT(out, check_code, symmetric, pivoting);
#else /* __TAHOE_MPI__ */
		/* constuctor */
		fLHS = new SPOOLESMatrixT(out, check_code, symmetric, pivoting);
#endif /* __TAHOE_MPI__ */
#else /* __SPOOLES__ */
		ExceptionT::GeneralFail(caller, "SPOOLES not installed");
#endif /* __SPOOLES__ */
	}
	else if (params.Name() == "SuperLU_matrix")
	{
		if (fFEManager.Size() == 1) /* serial */
		{
#ifdef __SUPERLU__
			/* global system properties */
			GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);

			bool symmetric;
			if (type == GlobalT::kDiagonal || type == GlobalT::kSymmetric)
				symmetric = true;
			else if (type == GlobalT::kNonSymmetric)
				symmetric = false;
			else
				ExceptionT::GeneralFail(caller, "unexpected system type: %d", type);

			/* construct */
			fLHS = new SuperLUMatrixT(out, check_code, symmetric);
#else /* no __SUPERLU__ */
			ExceptionT::GeneralFail(caller, "SuperLU not installed");
#endif /* __SUPERLU__*/
		}
		else /* parallel */
		{
#ifdef __SUPERLU_DIST__
			/* construct */
			fLHS = new SuperLU_DISTMatrixT(out, check_code, fFEManager.Communicator());
#else /* no __SUPERLU_DIST__ */
			ExceptionT::GeneralFail(caller, "SuperLU_DIST not installed");
#endif /* __SUPERLU_DIST__*/
		}
	}
	else if (params.Name() == "PSPASES_matrix")
	{
#ifdef __PSPASES__
		/* construct */
		fLHS = new PSPASESMatrixT(out, check_code, fFEManager.Communicator());
#else
		ExceptionT::GeneralFail(caller, " PSPASES solver not installed");
#endif /* __PSPASES__ */
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized matrix type \"%s\"", params.Name().Pointer());
		
#if 0
	switch (matrix_type)
	{
		case kAztec:
		{
#ifdef __AZTEC__
			/* construct */
			fLHS = new AztecMatrixT(in, out, check_code, fFEManager.Communicator());
#else
			ExceptionT::GeneralFail(caller, "Aztec solver not installed: %d", fMatrixType);
#endif /* __AZTEC__ */
			break;
		}
		default:
			ExceptionT::GeneralFail(caller, "unknown matrix type: %d", matrix_type);
	}	
#endif

	if (!fLHS) ExceptionT::OutOfMemory(caller);
}

/* call for equation renumbering */
bool SolverT::RenumberEquations(void)
{
	if (!fLHS) ExceptionT::GeneralFail("SolverT::RenumberEquations");
	return fLHS->RenumberEquations();
}