/* $Id: EVPFDBaseT.cpp,v 1.2.2.4 2001-07-02 21:54:34 paklein Exp $ */
/*
  File: EVPFDBaseT.cpp
*/

#include "EVPFDBaseT.h"
#include "NLCSolver.h"
#include "NLCSolver_LS.h"
#include "Utils.h"

#include "FEManagerT.h"
#include "FiniteStrainT.h"
#include "StringT.h"

/* initialization flag value */
const int kIsInit = 1;

/* spatial dimensions of the problem */
const int kNSD = 3;

EVPFDBaseT::EVPFDBaseT(ifstreamT& in, const FiniteStrainT& element) :
  FDHookeanMatT(in, element),
  IsotropicT  (in),
  //fdt         (element.FEManager().TimeStep()),
  ftime       (element.FEManager().Time()),
  fStatus     (element.RunState()),
  fLocDisp    (element.Displacements()),
  fLocLastDisp(element.LastDisplacements()),
  fKineticEqn (NULL),
  fSolver     (NULL),
  fSolverPtr  (new SolverWrapperEVPBase(*this)),
  fFtot       (kNSD),
  fs_ij       (kNSD),
  fc_ijkl     (dSymMatrixT::NumValues(kNSD))
{
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

bool EVPFDBaseT::NeedsInitialization() const { return false; }
void EVPFDBaseT::Initialize()
{
  // set slip system hardening law
  SetKineticEquation();

  // set nonlinear constitutive solver 
  SetConstitutiveSolver();

  // allocate space for all elements
  AllocateElements();

  // initialize selected variables in all elements
  InitializeVariables();
}

void EVPFDBaseT::Print(ostream& out) const
{
  // inherited
  FDHookeanMatT::Print(out);
  IsotropicT::Print(out);

  // print material data
  out << " Material data:\n";
  out << "    Kinetics of deformation\n";
  out << "       Kinetic eqn model code. . . . . . . . . . = " << fKinEqnCode << "\n";
  fKineticEqn->Print(out);

  // nonlinear constitutive solver data
  out << "    NLC solver:\n";
  out << "       Solver method . . . . . . . . . . . . . . = " << fSolverCode << "\n";
  fSolver->Print(out);
}

/* set (material) tangent modulus */
void EVPFDBaseT::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

void EVPFDBaseT::PrintName(ostream& out) const
{
  // inherited
  FDHookeanMatT::PrintName(out);

  // print model name
  out << "    Elasto-Visco-Plastic constitutive model\n";
}

int EVPFDBaseT::GetNumberOfEqns()
{
  // must be defined in derived classes
  throwRunTimeError("EVPFDBaseT::GetNumberOfEqns: called!!!");
  return 0;
}

void EVPFDBaseT::AllocateElements()
{
  // determine storage size
  int i_size = NumIP();
  int d_size = NumVariablesPerElement();

  // allocate space for all elements
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);
      
      // construct element
      element.Allocate(i_size, d_size);
      
      // initialize values
      element.IntegerData() = kIsInit;
      element.DoubleData()  = 0.0;
    }
}

// compute 3D deformation gradient
void EVPFDBaseT::Compute_Ftot_3D(dMatrixT& F_3D) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D =  F();
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F());    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw eGeneralFail;
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
