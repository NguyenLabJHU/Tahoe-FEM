/* $Id: EVPFDBaseT.cpp,v 1.12 2004-07-15 08:29:14 paklein Exp $ */
#include "EVPFDBaseT.h"
#include "NLCSolver.h"
#include "NLCSolver_LS.h"
#include "Utils.h"

#include "FSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* initialization flag value */
const int kIsInit = 1;

/* spatial dimensions of the problem */
const int kNSD = 3;

EVPFDBaseT::EVPFDBaseT(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("EVPFDBase"),
//  FDHookeanMatT(in, support),
  IsotropicT  (in),
  //fdt         (element.FEManager().TimeStep()),
  //ftime       (element.ElementSupport().Time()),
  //fStatus     (element.RunState()),
  fLocDisp    (support.LocalArray(LocalArrayT::kDisp)),
  fLocLastDisp(support.LocalArray(LocalArrayT::kLastDisp)),
  fKineticEqn (NULL),
  fSolver     (NULL),
  fSolverPtr  (new SolverWrapperEVPBase(*this)),
  fFtot       (kNSD),
  fs_ij       (kNSD),
  fc_ijkl     (dSymMatrixT::NumValues(kNSD))
{
ExceptionT::GeneralFail("EVPFDBaseT::EVPFDBaseT", "out of date");

  // input file
  StringT filename;
  in >> filename;
  
  // generate relative path in native format
  filename.ToNativePathName();
  StringT path;
  path.FilePath(in.filename());
  filename.Prepend(path);
  
  OpenExternal(fInput, filename, "EVPFDBaseT data");
  if (in.skip_comments())
    fInput.set_marker(in.comment_marker());

  // Lame's constants
  fmu     = Mu();
  flambda = Lambda();
  fbulk   = flambda + 2./3.*fmu;
}

EVPFDBaseT::~EVPFDBaseT()
{
  delete fKineticEqn;
  delete fSolver;
}


void EVPFDBaseT::Initialize()
{
  // set slip system hardening law
  SetKineticEquation();

  // set nonlinear constitutive solver 
  SetConstitutiveSolver();

  // allocate space for all elements
  //AllocateElements();

  // initialize selected variables in all elements
  //InitializeVariables();
}

/* apply initialize current time step */
void EVPFDBaseT::InitStep(void)
{
	/* inherited */
	FDHookeanMatT::InitStep();
	
	/* cache the time step */
	fdt = MaterialSupport().TimeStep();
}

bool EVPFDBaseT::NeedLastDisp() const { return true; }

/* information about subordinate parameter lists */
void EVPFDBaseT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FDHookeanMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EVPFDBaseT::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = FDHookeanMatT::NewSub(name);
	if (params)
		return params;
	else
		return IsotropicT::NewSub(name);
}

/* accept parameter list */
void EVPFDBaseT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FDHookeanMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
}

/* set (material) tangent modulus */
void EVPFDBaseT::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

int EVPFDBaseT::GetNumberOfEqns()
{
  // must be defined in derived classes
  throwRunTimeError("EVPFDBaseT::GetNumberOfEqns: called!!!");
  return 0;
}

bool EVPFDBaseT::NeedsPointInitialization() const { return true; }
void EVPFDBaseT::PointInitialize()
{
	// only execute during the first ip
	if (CurrIP() == 0) {

		// determine storage size
		int i_size = NumIP();
		int d_size = NumVariablesPerElement();

		// current element
		ElementCardT& element = CurrentElement();	

		// construct element
		element.Dimension(i_size, d_size);
		element.IntegerData() = kIsInit;
		element.DoubleData()  = 0.0;		

		// initialize values
		InitializeVariables(element);
    }
}

// compute 3D deformation gradient
void EVPFDBaseT::Compute_Ftot_3D(dMatrixT& F_3D) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D = F();
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F());
		F_3D(2, 2) = 1.0;
	}
	else 
		throw ExceptionT::kGeneralFail;
}

void EVPFDBaseT::Compute_Ftot_last_3D(dMatrixT& F_3D) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D = F_total_last();
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F_total_last());
		F_3D(2, 2) = 1.0;
	}
	else 
		throw ExceptionT::kGeneralFail;
}

void EVPFDBaseT::SetConstitutiveSolver()
{
  // input solver code
  fInput >> fSolverCode;

  // get size of nonlinear system of equations
  int numeqs = GetNumberOfEqns();

  // select constitutive solver
  switch(fSolverCode)
    {
      // Newton's method + line search
    case NLCSolver::kNLCSolver_LS: 
      fSolver = new NLCSolver_LS(numeqs);
      break;

    default:
      throwRunTimeError("EVPFDBaseT::SetConstitutiveSolver: Bad fSolverCode");
    }
  if (!fSolver) throwMemoryError("EVPFDBaseT::SetConstitutiveSolver");

  // modify some default values of solver
  int maxiter;
  fInput >> maxiter;
  fSolver->SetMaxIterations(maxiter);

  double functol;
  fInput >> functol;
  fSolver->SetFuncTol(functol);

  double gradtol;
  fInput >> gradtol;
  fSolver->SetGradTol(gradtol);
}
