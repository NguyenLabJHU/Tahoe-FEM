/* $Id: SolverT.cpp,v 1.39 2007-01-14 22:43:25 paklein Exp $ */
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
#include "SuperLUMatrixT.h"
#include "SuperLU_MTMatrixT.h"
#include "SuperLU_DISTMatrixT.h"
#include "SPOOLESMatrixT.h"
#include "PSPASESMatrixT.h"

#ifdef __AZTEC__
#include "AztecParamsT.h"
#include "AztecMatrixT.h"
#endif

#ifdef __TAHOE_MPI__
#include "SPOOLESMatrixT_mpi.h"
#endif

#ifdef __SPOOLES_MT__
#include "SPOOLESMatrixT_MT.h"
#endif

#ifdef __TRILINOS__
/* Trilinos-Anasazi headers */
#define HAVE_CONFIG_H
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include "TrilinosAztecT.h"
#include "NodeManagerT.h"
#include "FieldT.h"
#include "nIntegratorT.h"
#include "eStaticIntegrator.h"
using namespace Anasazi;

/* translate integer to Anasazi message type */
static Anasazi::MsgType int2MsgType(int i) {
	switch (i) {
		case Anasazi::Warnings:
			return Anasazi::Warnings;
		case Anasazi::IterationDetails:
			return Anasazi::IterationDetails;
		case Anasazi::OrthoDetails:
			return Anasazi::OrthoDetails;
		case Anasazi::FinalSummary:
			return Anasazi::FinalSummary;
		case Anasazi::TimingDetails:
			return Anasazi::TimingDetails;
		case Anasazi::StatusTestDetails:
			return Anasazi::StatusTestDetails;
		case Anasazi::Debug:
			return Anasazi::Debug;
	}
	return Anasazi::Errors;
};
#endif /* __TRILINOS__ */

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
	fPerturbation(0.0),
	fEigenSolverParameters(NULL)
{
	/* console */
	iSetName("solver");
	iAddVariable("print_equation_numbers", fPrintEquationNumbers);
}

/* destructor */
SolverT::~SolverT(void) {
	delete fLHS;
	delete fEigenSolverParameters;
}

/* configure the global equation system */
void SolverT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{	
	try {
// DEBUG
//cout << "\n" << "Total # of equations: " << tot_num_eq << "\n" << "Local # of equations: " << loc_num_eq << "\n" << "First eq on this proc: " << start_eq << endl;
		
		/* allocate rhs vector */
		fRHS.Dimension(loc_num_eq);
		fRHS = 0.0;
		
		/* set global equation matrix type */
		fLHS->Initialize(tot_num_eq, loc_num_eq, start_eq);
		
		/* write information */
		if (fFEManager.Logging() != GlobalT::kSilent)
			fLHS->Info(fFEManager.Output());
	
		/* output global equation number for each DOF */
		if (fPrintEquationNumbers) fFEManager.WriteEquationNumbers(fGroup);
	}	

	catch (ExceptionT::CodeT error_code) {
		ExceptionT::Throw(error_code, "SolverT::Initialize");
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
	const char caller[] = "SolverT::CloseStep";

#ifdef __TRILINOS__
	if (!fEigenSolverParameters) return; /* no eigenmode solver defined */

	/* temporarily replace LHS */
	EpetraCRSMatrixT epetra(fFEManager.Output(), GlobalMatrixT::kNoCheck, fLHS->Communicator());
	GlobalMatrixT* p_lhs = fLHS;
	fLHS = &epetra;
	
	/* set storage */
	fFEManager.SendEqnsToSolver(fGroup);	
	epetra.Initialize(p_lhs->NumTotEquations(), p_lhs->NumEquations(), p_lhs->StartEquation());
	
	/* get time integrator */
	const NodeManagerT* nodes = fFEManager.NodeManager();
	ArrayT<FieldT*> fields;
	nodes->CollectFields(fGroup, fields);
	
	/* calculate K only */
	for (int i = 0; i < fields.Length(); i++) {
		eIntegratorT& e_int = const_cast<eIntegratorT&>(fields[i]->nIntegrator().eIntegrator());
		eStaticIntegrator* e_static = TB_DYNAMIC_CAST(eStaticIntegrator*, &e_int);
		if (! e_static) ExceptionT::GeneralFail(caller, "Could not cast integrator to eStaticIntegrator for field \"%s\"", fields[i]->FieldName().Pointer());
		e_static->SetLHSMode(eStaticIntegrator::kFormKOnly);
	}

	/* calc K and store */
	fLHS->Clear();
	fLHS_lock = kOpen;
	fFEManager.FormLHS(Group(), fLHS->MatrixType());
	fLHS_lock = kLocked;
	Epetra_CrsMatrix* K = epetra.Translate();
	
	/* calculate M only */
	for (int i = 0; i < fields.Length(); i++) {
		eIntegratorT& e_int = const_cast<eIntegratorT&>(fields[i]->nIntegrator().eIntegrator());
		eStaticIntegrator* e_static = TB_DYNAMIC_CAST(eStaticIntegrator*, &e_int);
		if (! e_static) ExceptionT::GeneralFail(caller, "Could not cast integrator to eStaticIntegrator for field \"%s\"", fields[i]->FieldName().Pointer());
		e_static->SetLHSMode(eStaticIntegrator::kFormMOnly);
	}
	
	/* calc M and store */
	fLHS->Clear();
	fLHS_lock = kOpen;
	fFEManager.FormLHS(Group(), fLHS->MatrixType());
	fLHS_lock = kLocked;
	Epetra_CrsMatrix* M = epetra.Translate();

	/* restore LHS */
	fLHS = p_lhs;
	
	/* restore integrators */
	for (int i = 0; i < fields.Length(); i++) {
		eIntegratorT& e_int = const_cast<eIntegratorT&>(fields[i]->nIntegrator().eIntegrator());
		eStaticIntegrator* e_static = TB_DYNAMIC_CAST(eStaticIntegrator*, &e_int);
		if (! e_static) ExceptionT::GeneralFail(caller, "Could not cast integrator to eStaticIntegrator for field \"%s\"", fields[i]->FieldName().Pointer());
		e_static->SetLHSMode(eStaticIntegrator::kNormal);
	}
	
	/* set up eigensolver */
	Teuchos::RefCountPtr<Epetra_CrsMatrix> rcp_K = Teuchos::rcp(K);
	Teuchos::RefCountPtr<Epetra_CrsMatrix> rcp_M = Teuchos::rcp(M);
	
	/* write matricies to files */
	bool output_matricies = fEigenSolverParameters->GetParameter("output_matricies");
	if (output_matricies) {
		ofstreamT out;
		out.precision(12);

		StringT root = fFEManager.InputFile();
		root.Root();

		StringT k_file = root;
		k_file.Append(".K.", fFEManager.StepNumber(), 2);
		k_file.Append(".dat");
		out.open(k_file);
		K->Print(out);
		out.close();
		
		StringT m_file = root;
		m_file.Append(".M.", fFEManager.StepNumber(), 2);
		m_file.Append(".dat");
		out.open(m_file);
		M->Print(out);
		out.close();
	}

	/* parameters common to all solvers */
	int nev = fEigenSolverParameters->GetParameter("max_eigenmodes");
	int block_size = fEigenSolverParameters->GetParameter("block_size");

	/* work space */
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

	// Create an Epetra_MultiVector for an initial vector to start the solver.
	// Note:  This needs to have the same number of columns as the blocksize.
	Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(*(epetra.Map()), block_size) );
	ivec->Random();

	// Create the eigenproblem.
	Teuchos::RefCountPtr<BasicEigenproblem<double, MV, OP> > MyProblem =
		Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(rcp_K, rcp_M, ivec) );

	/* global system properties */
	GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);
	if (type == GlobalT::kSymmetric)
		MyProblem->setHermitian(true);
	else
		MyProblem->setHermitian(false);

	// Set the number of eigenvalues requested
	MyProblem->setNEV(nev);

	// Inform the eigenproblem that you are finishing passing it information
	bool boolret = MyProblem->setProblem();
	if (boolret != true) {
		ExceptionT::GeneralFail(caller, "Anasazi::BasicEigenproblem::setProblem() returned an error");
	}

	/* which [in] The eigenvalues of interest for this eigenproblem.
	* "LM" - Largest Magnitude [ default ]
	* "SM" - Smallest Magnitude
	* "LR" - Largest Real
	* "SR" - Smallest Real
	* "LI" - Largest Imaginary
	* "SI" - Smallest Imaginary */
	std::string which("SM"); /* get lowest modes */

	/* Create parameter list to pass into the solver manager */
	Teuchos::ParameterList MyPL;
	MyPL.set("Which", which);
	MyPL.set( "Block Size", block_size);
	Anasazi::SolverManager<double, MV, OP>* solver_man = NULL;
	
	/* construct solver */
	if (fEigenSolverParameters->Name() == "LOBPCG_solver") {

		/* resolve parameters */
		const int max_iter = fEigenSolverParameters->GetParameter("max_iter");
		const double tol   = fEigenSolverParameters->GetParameter("tolerance");
		Anasazi::MsgType verbosity = int2MsgType(fEigenSolverParameters->GetParameter("verbosity"));

		/* LOBPCG parameters */
		MyPL.set("Maximum Iterations", max_iter);
		MyPL.set("Convergence Tolerance", tol);
		MyPL.set("Verbosity", verbosity);

		/* Create the solver manager */
		solver_man = new SimpleLOBPCGSolMgr<double, MV, OP>(MyProblem, MyPL);
	
	} else if (fEigenSolverParameters->Name() == "Block_Davidson_solver") {

		/* resolve parameters */
		const int max_restart = fEigenSolverParameters->GetParameter("max_restart");
		const int num_blocks = fEigenSolverParameters->GetParameter("num_blocks");
		const double tol   = fEigenSolverParameters->GetParameter("tolerance");
		Anasazi::MsgType verbosity = int2MsgType(fEigenSolverParameters->GetParameter("verbosity"));

		/* Block Davidson parameters */
		MyPL.set("Num Blocks", num_blocks);
		MyPL.set("Maximum Restarts", max_restart);
		MyPL.set("Convergence Tolerance", tol);
		MyPL.set("Verbosity", verbosity);

		/* Create the solver manager */
		solver_man = new Anasazi::BlockDavidsonSolMgr<double, MV, OP>(MyProblem, MyPL);

	} else
		ExceptionT::GeneralFail(caller, "unrecognized solver name [%s]", fEigenSolverParameters->Name().Pointer());

	/* Solve the problem */
	try {
		ReturnType returnCode = solver_man->solve();
		if (returnCode != Anasazi::Converged) {
			cout << caller << ": eigensystem solve did not converge" << endl;
		}

		// Get the eigenvalues and eigenvectors from the eigenproblem
		Eigensolution<double,MV> sol = MyProblem->getSolution();
		std::vector<Value<double> > evals = sol.Evals;
		Teuchos::RefCountPtr<MV> evecs = sol.Evecs;

		for (int i = 0; i <  sol.numVecs; i++) {
			double re = evals[i].realpart;
			double im = evals[i].imagpart;
			cout << i << ": " << re << ", " << im << '\n';
		}
		cout.flush();
	}
	catch (std::exception e) {
		cout << caller << ": caught exception from Trilinos" << endl;
	}
	
	// set up output
	
	// loop over modes and write modes and write e-vals to .out
	
	/* clean up */
	delete solver_man;
//	delete K;
//	delete M;
#endif
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
		ExceptionT::GeneralFail("SolverT::EquationNumberScope", "invalid LHS");
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
	check_code.AddEnumeration("all_pivots", GlobalMatrixT::kAllPivots);
	check_code.AddEnumeration("print_LHS", GlobalMatrixT::kPrintLHS);
	check_code.AddEnumeration("print_RHS", GlobalMatrixT::kPrintRHS);
	check_code.AddEnumeration("print_solution", GlobalMatrixT::kPrintSolution);
	check_code.AddEnumeration("check_LHS", GlobalMatrixT::kCheckLHS);
#ifdef __TRILINOS__
	check_code.AddEnumeration("eigenmodes", GlobalMatrixT::kEigenmodes);
#endif
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

#ifdef __TRILINOS__
	/* eigenvalue solver options */
	sub_list.AddSub("eigensolver_choice", ParameterListT::ZeroOrOnce, true);
#endif
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SolverT::NewSub(const StringT& name) const
{
	if (name == "matrix_type_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		choice->AddSub(ParameterContainerT("profile_matrix"));
		choice->AddSub(ParameterContainerT("diagonal_matrix"));
		choice->AddSub(ParameterContainerT("full_matrix"));

#ifdef __AZTEC__
		choice->AddSub("Aztec_matrix");
#endif

#ifdef __TRILINOS__
		ParameterContainerT Trilinos_Aztec("Trilinos_Aztec_matrix");

		/* options */
		Trilinos_Aztec.AddParameter(ParameterT::Integer, "dummy");

		choice->AddSub(Trilinos_Aztec);
#endif

#ifdef __SPOOLES__
		ParameterContainerT SPOOLES("SPOOLES_matrix");
		ParameterT message_level(ParameterT::Enumeration, "message_level");
		message_level.AddEnumeration("silent", 0);
		message_level.AddEnumeration("timing", 1);
		message_level.AddEnumeration("verbose", 99);
		message_level.SetDefault(0);
		SPOOLES.AddParameter(message_level);
		ParameterT enable_pivoting(ParameterT::Boolean, "enable_pivoting");
		enable_pivoting.SetDefault(true);
		SPOOLES.AddParameter(enable_pivoting);
		ParameterT always_symmetric(ParameterT::Boolean, "always_symmetric");
		always_symmetric.SetDefault(false);
		SPOOLES.AddParameter(always_symmetric);
		choice->AddSub(SPOOLES);

#ifdef __SPOOLES_MT__
		ParameterContainerT SPOOLES_MT(SPOOLES);
		SPOOLES_MT.SetName("SPOOLES_MT_matrix");		
		ParameterT num_threads(ParameterT::Integer, "num_threads");
		num_threads.AddLimit(2, LimitT::LowerInclusive);
		SPOOLES_MT.AddParameter(num_threads);
		choice->AddSub(SPOOLES_MT);
#endif /* __SPOOLES_MT__ */
#endif /* __SPOOLES__ */

#ifdef __PSPASES__
		choice->AddSub(ParameterContainerT("PSPASES_matrix"));
#endif	

#if defined(__SUPERLU__) || defined(__SUPERLU_DIST__)
		/* output timing statistics */
		ParameterT print_stat(ParameterT::Boolean, "print_stat");
		print_stat.SetDefault(false);
		
		/* solution refinement */
		ParameterT refinement(ParameterT::Enumeration, "refinement");
		refinement.AddEnumeration("NOREFINE", NOREFINE);
		refinement.AddEnumeration(  "SINGLE", SINGLE);
		refinement.AddEnumeration(  "DOUBLE", DOUBLE);
		refinement.AddEnumeration(   "EXTRA", EXTRA);
		refinement.SetDefault(NOREFINE);
#endif /* __SUPERLU__ or __SUPERLU_DIST__ */

#ifdef __SUPERLU__
		ParameterContainerT SuperLU("SuperLU_matrix");

		/* options */
		SuperLU.AddParameter(print_stat);
		SuperLU.AddParameter(refinement);		

		choice->AddSub(SuperLU);
#endif /* __SUPERLU__ */

#ifdef __SUPERLU_DIST__
		ParameterContainerT SuperLU_DIST("SuperLU_DIST_matrix");

		/* options */
		SuperLU_DIST.AddParameter(print_stat);
		SuperLU_DIST.AddParameter(refinement);		

		choice->AddSub(SuperLU_DIST);
#endif /* __SUPERLU_DIST__ */

#ifdef __SUPERLU_MT__
		ParameterContainerT SuperLU_MT("SuperLU_MT_matrix");
		ParameterT slu_num_threads(ParameterT::Integer, "num_threads");
		slu_num_threads.AddLimit(1, LimitT::LowerInclusive);
		SuperLU_MT.AddParameter(slu_num_threads);
		choice->AddSub(SuperLU_MT);
#endif /* __SUPERLU_MT__ */

		return choice;
	}
#ifdef __AZTEC__
	else if (name == "Aztec_matrix")
		return new AztecParamsT;
#endif

#ifdef __TRILINOS__
	else if (name == "eigensolver_choice") {
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		/* number of eigenmodes */
		ParameterT num_modes(ParameterT::Integer, "max_eigenmodes");
		num_modes.AddLimit(0, LimitT::LowerInclusive);
		num_modes.SetDefault(10);

		/* block size */
		ParameterT block_size(ParameterT::Integer, "block_size");
		block_size.AddLimit(0, LimitT::LowerInclusive);
		block_size.SetDefault(5);

		/* convergence tolerance */
		ParameterT tol(ParameterT::Double, "tolerance");
		tol.SetDefault(1.0e-10);

		/* write matricies to files */
		ParameterT output(ParameterT::Boolean, "output_matricies");
		output.SetDefault(false);

		/* solver verbosity */
		ParameterT verbosity(ParameterT::Enumeration, "verbosity");
		verbosity.AddEnumeration(           "Errors", Anasazi::Errors);
		verbosity.AddEnumeration(         "Warnings", Anasazi::Warnings);
		verbosity.AddEnumeration( "IterationDetails", Anasazi::IterationDetails);
		verbosity.AddEnumeration(     "OrthoDetails", Anasazi::OrthoDetails);
		verbosity.AddEnumeration(     "FinalSummary", Anasazi::FinalSummary);
		verbosity.AddEnumeration(    "TimingDetails", Anasazi::TimingDetails);
		verbosity.AddEnumeration("StatusTestDetails", Anasazi::StatusTestDetails);
		verbosity.AddEnumeration(            "Debug", Anasazi::Debug);
		verbosity.SetDefault(Anasazi::Errors);
// from AnasaziTypes.hpp
#if 0
  enum MsgType 
  {
    Errors = 0,                 /*!< Errors [ always printed ] */
    Warnings = 0x1,             /*!< Internal warnings */
    IterationDetails = 0x2,     /*!< Approximate eigenvalues, errors */
    OrthoDetails = 0x4,         /*!< Orthogonalization/orthonormalization details */
    FinalSummary = 0x8,         /*!< Final computational summary */
    TimingDetails = 0x10,       /*!< Timing details */
    StatusTestDetails = 0x20,   /*!< Status test details */
    Debug = 0x40                /*!< Debugging information */
  };
#endif

		/* LOBPCG Method */
		ParameterContainerT LOBPCG("LOBPCG_solver");
		LOBPCG.AddParameter(num_modes);
		LOBPCG.AddParameter(block_size);
		LOBPCG.AddParameter(tol);
		LOBPCG.AddParameter(output);
		LOBPCG.AddParameter(verbosity);

		ParameterT max_iter(ParameterT::Integer, "max_iter");
		max_iter.AddLimit(0, LimitT::LowerInclusive);
		max_iter.SetDefault(500);
		LOBPCG.AddParameter(max_iter);
		
		choice->AddSub(LOBPCG);

		/* Block Davidson Method */
		ParameterContainerT BlockDavidson("Block_Davidson_solver");
		BlockDavidson.AddParameter(num_modes);
		BlockDavidson.AddParameter(block_size);
		BlockDavidson.AddParameter(tol);
		BlockDavidson.AddParameter(output);		
		BlockDavidson.AddParameter(verbosity);

		ParameterT max_restart(ParameterT::Integer, "max_restart");
		max_restart.AddLimit(0, LimitT::LowerInclusive);
		max_restart.SetDefault(500);
		BlockDavidson.AddParameter(max_restart);

		ParameterT num_blocks(ParameterT::Integer, "num_blocks");
		num_blocks.AddLimit(0, LimitT::LowerInclusive);
		num_blocks.SetDefault(8);
		BlockDavidson.AddParameter(num_blocks);

		choice->AddSub(BlockDavidson);
		
		return choice;
	}
#endif

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

#ifdef __TRILINOS__
	/* look for eigenmodes solver */
	const ParameterListT* eig_solve = list.ListChoice(*this, "eigensolver_choice");
	if (eig_solve) {
		delete fEigenSolverParameters;
		fEigenSolverParameters = new ParameterListT(*eig_solve);
	}
#endif
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
	//approx_LHS->SetPrintTag("approximate_LHS.");

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

	out << "\nreference LHS (even suffix):\n";
	ref_LHS.PrintLHS(true);

	out << "\ntest LHS (odd suffix):\n";
	test_LHS.PrintLHS(true);
	
	out.flush();
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
	ofstreamT& out = fFEManager.Output();
	
	/* MP support */
	const CommunicatorT& comm = fFEManager.Communicator();

	/* resolve matrix type */
	if (params.Name() == "diagonal_matrix")
	{
		DiagonalMatrixT* diag = new DiagonalMatrixT(out, check_code, DiagonalMatrixT::kNoAssembly, comm);
		diag->SetAssemblyMode(DiagonalMatrixT::kDiagOnly);
		fLHS = diag;
	}
	else if (params.Name() == "profile_matrix")
	{
		/* global system properties */
		GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);
		
		if (type == GlobalT::kNonSymmetric)
			fLHS = new CCNSMatrixT(out, check_code, comm);
		else if (type == GlobalT::kSymmetric)
			fLHS = new CCSMatrixT(out, check_code, comm);
		else
			ExceptionT::GeneralFail(caller, "system type %d not compatible with matrix %d", type, kProfileSolver);
	}
	else if (params.Name() == "full_matrix")
		fLHS = new FullMatrixT(out, check_code, comm);
	else if (params.Name() == "SPOOLES_matrix" || params.Name() == "SPOOLES_MT_matrix")
	{
#ifdef __SPOOLES__
		/* global system properties */
		GlobalT::SystemTypeT type = fFEManager.GlobalSystemType(fGroup);

		/* solver options */
		int message_level = params.GetParameter("message_level");
		bool pivoting = params.GetParameter("enable_pivoting");
		bool always_symmetric = params.GetParameter("always_symmetric");
		bool symmetric;
		if (always_symmetric)
			symmetric = true;
		else if (type == GlobalT::kDiagonal || type == GlobalT::kSymmetric)
			symmetric = true;
		else if (type == GlobalT::kNonSymmetric)
			symmetric = false;
		else
			ExceptionT::GeneralFail(caller, "unexpected system type: %d", type);
		
		/* check */
		if (!symmetric && !pivoting)
			ExceptionT::GeneralFail(caller, "pivoting required with non-symmetric matricies");
		// NOTE: SPOOLES v2.2 does not seem to solve non-symmetric
		//      systems correctly in parallel if pivoting is disabled

		/* multi-threaded SPOOLES */
		if (params.Name() == "SPOOLES_MT_matrix")
		{
			/* number of solver threads */
			int num_threads = params.GetParameter("num_threads");
		
#ifdef __SPOOLES_MT__
				fLHS = new SPOOLESMatrixT_MT(out, check_code, symmetric, pivoting, message_level, num_threads, comm);
#else /* __SPOOLES_MT__ */
				ExceptionT::GeneralFail(caller, "SPOOLES MPI not installed");
#endif /* __SPOOLES_MT__ */
		}
		else {
#ifdef __TAHOE_MPI__
			/* constuctor */
			if (fFEManager.Size() > 1)
			{
#ifdef __SPOOLES_MPI__
				fLHS = new SPOOLESMatrixT_mpi(out, check_code, symmetric, pivoting, message_level, comm);
#else /* __SPOOLES_MPI__ */
				ExceptionT::GeneralFail(caller, "SPOOLES MPI not installed");
#endif /* __SPOOLES_MPI__ */
			}
			else /* single processor with MPI-enabled code */
				fLHS = new SPOOLESMatrixT(out, check_code, symmetric, pivoting, message_level, comm);
#else /* __TAHOE_MPI__ */
			/* constuctor */
			fLHS = new SPOOLESMatrixT(out, check_code, symmetric, pivoting, message_level, comm);
#endif /* __TAHOE_MPI__ */
		}
#else /* __SPOOLES__ */
		ExceptionT::GeneralFail(caller, "SPOOLES not installed");
#endif /* __SPOOLES__ */
	}
	else if (params.Name() == "SuperLU_matrix")
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
		bool print_stat = params.GetParameter("print_stat");
		int i_refine  = params.GetParameter("refinement");		
		IterRefine_t refine = SuperLUMatrixT::int2IterRefine_t(i_refine);

		fLHS = new SuperLUMatrixT(out, check_code, symmetric, print_stat, refine, comm);
#else /* no __SUPERLU__ */
		ExceptionT::GeneralFail(caller, "SuperLU not installed");
#endif /* __SUPERLU__*/
	}
	else if (params.Name() == "SuperLU_DIST_matrix")
	{
#ifdef __SUPERLU_DIST__
		/* construct */
		bool print_stat = params.GetParameter("print_stat");
		int i_refine  = params.GetParameter("refinement");
		IterRefine_t refine = SuperLU_DISTMatrixT::int2IterRefine_t(i_refine);

		fLHS = new SuperLU_DISTMatrixT(out, check_code, print_stat, refine, comm);
#else /* no __SUPERLU_DIST__ */
		ExceptionT::GeneralFail(caller, "SuperLU_DIST not installed");
#endif /* __SUPERLU_DIST__ */
	}
	else if (params.Name() == "SuperLU_MT_matrix")
	{
#ifdef __SUPERLU_MT__
		int num_threads = params.GetParameter("num_threads");
		fLHS = new SuperLU_MTMatrixT(out, check_code, num_threads, comm);
#else
		ExceptionT::GeneralFail(caller, "SuperLU_MT not installed");
#endif	
	}
	else if (params.Name() == "PSPASES_matrix")
	{
#ifdef __PSPASES__
		/* construct */
		fLHS = new PSPASESMatrixT(out, check_code, comm);
#else
		ExceptionT::GeneralFail(caller, " PSPASES solver not installed");
#endif /* __PSPASES__ */
	}
	else if (params.Name() == "Aztec_matrix")
	{
#ifdef __AZTEC__
			/* construct */
			fLHS = new AztecMatrixT(out, check_code, comm, params);
#else
			ExceptionT::GeneralFail(caller, "Aztec solver not installed: %d", fMatrixType);
#endif /* __AZTEC__ */
	}
	else if (params.Name() == "Trilinos_Aztec_matrix")
	{
#ifdef __TRILINOS__
			/* construct */
			fLHS = new TrilinosAztecT(out, check_code, comm);
#else
			ExceptionT::GeneralFail(caller, "Trilinos not installed: %d", fMatrixType);
#endif /* __TRILINOS__ */
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized matrix type \"%s\"", params.Name().Pointer());

	if (!fLHS) ExceptionT::OutOfMemory(caller);
}

/* call for equation renumbering */
bool SolverT::RenumberEquations(void)
{
	if (!fLHS) ExceptionT::GeneralFail("SolverT::RenumberEquations");
	return fLHS->RenumberEquations();
}
