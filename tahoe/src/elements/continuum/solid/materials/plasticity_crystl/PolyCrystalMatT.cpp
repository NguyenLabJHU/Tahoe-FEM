/*
  File: PolyCrystalMatT.cpp
*/

#include "PolyCrystalMatT.h"
#include "CrystalElasticity.h"
#include "SlipGeometry.h"
#include "LatticeOrient.h"
#include "NLCSolver.h"
#include "NLCSolver_LS.h"
#include "SlipKinetics.h"
#include "SlipHardening.h"
#include "Utils.h"

#include "FEManagerT.h"
#include "FiniteStrainT.h"
#include "StringT.h"

/* number of elastic material properties : isotropic and cubic */
const int kNumMatProp = 3;

/* number of slip systems based on crystal structure */
const int kSlipFCC = 12;
const int kSlipBCC = 12;
const int kSlipHCP = 12;

/* initialization flag value */
const int kIsInit = 1;

/* spatial dimensions of the problem */
const int kNSD = 3;

PolyCrystalMatT::PolyCrystalMatT(ifstreamT& in, const FiniteStrainT& element) :
  FDHookeanMatT(in, element),
  //fdt           (element.FEManager().TimeStep()),
  ftime         (element.FEManager().Time()),
  fStatus       (element.RunState()),
  fLocLastDisp  (element.LastDisplacements()),
  fLocDisp      (element.Displacements()),
  fSlipGeometry (NULL),
  fLatticeOrient(NULL),
  fElasticity   (NULL),
  fKinetics     (NULL),
  fHardening    (NULL),
  fSolver       (NULL),
  fSolverPtr    (new SolverWrapperPoly(*this)),
  fMatProp      (kNumMatProp),
  fFtot_n       (kNSD,kNSD),
  fFtot         (kNSD,kNSD),
  fsavg_ij      (kNSD),
  fcavg_ijkl    (dSymMatrixT::NumValues(kNSD))
{
  // input file
  StringT filename;
  in >> filename;
  
  // generate relative path in native format
  filename.ToNativePathName();
  StringT path;
  path.FilePath(in.filename());
  filename.Prepend(path);

  OpenExternal(fInput, filename, "PolyCrystalMatT data");
  if (in.skip_comments())
    fInput.set_marker(in.comment_marker());

  // read number of crystals per integration point 
  fInput >> fNumGrain;

  // set slip system quantities
  SetSlipSystems();

  // read crystal orientations
  SetLatticeOrientation();

  // set crystal elasticity type
  SetCrystalElasticity();

  // set nonlinear constitutive solver 
  SetConstitutiveSolver();

  // read data for state iteration
  fInput >> fMaxIterState;
  fInput >> fTolerState;
}

PolyCrystalMatT::~PolyCrystalMatT()
{
  delete fSlipGeometry;
  delete fLatticeOrient;
  delete fElasticity;
  delete fKinetics;
  delete fHardening;
  delete fSolver;
}

bool PolyCrystalMatT::NeedsInitialization() const { return false; }
void PolyCrystalMatT::Initialize()
{
  // set slip system hardening law (set before slip kinetics!!!)
  SetSlipHardening();

  // set kinetics of slip
  SetSlipKinetics();

  // allocate space for all elements
  AllocateElements();

  // initialize state variables in all elements
  InitializeCrystalVariables();
}

bool PolyCrystalMatT::NeedLastDisp() const { return true; }

void PolyCrystalMatT::Print(ostream& out) const
{
  // inherited
  FDHookeanMatT::Print(out);

  // print input values
  out << " Polycrystal data:\n";
  out << "    Number of grains   . . . . . . . . . . . . . = " << fNumGrain    << "\n";
  out << "    Crystal structure  . . . . . . . . . . . . . = " << fCrystalType << "\n";
  out << "    Lattice Orientation code . . . . . . . . . . = " << fODFCode << "\n";
  out << "    Incrs to output ODF. . . . . . . . . . . . . = " << fODFOutInc << "\n";
 
  // elasticity constants
  out << "    Elasticity type. . . . . . . . . . . . . . . = " << fElastCode << "\n";
  fElasticity->Print(out);

  // control data for gamma solver
  out << "    NLC solver for dgamma\n";
  out << "       Solver method . . . . . . . . . . . . . . = " << fSolverCode << "\n";
  fSolver->Print(out);

  // iterations on state
  out << "    Convergence control for state\n";
  out << "       Max# iterations . . . . . . . . . . . . . = " << fMaxIterState << "\n";
  out << "       Tolerance convergence . . . . . . . . . . = " << fTolerState   << "\n";
}

/* set (material) tangent modulus */
void PolyCrystalMatT::SetModulus(dMatrixT& modulus)
{
	fElasticity->ComputeModuli(modulus);
}

void PolyCrystalMatT::PrintName(ostream& out) const
{
  // inherited
  FDHookeanMatT::PrintName(out);

  // print model name
  out << "    PolyCrystal constitutive model\n";
}

void PolyCrystalMatT::AllocateElements()
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

void PolyCrystalMatT::SetSlipSystems()
{
  // read crystal structure
  fInput >> fCrystalType;

  // construct slip system geometry
  switch(fCrystalType)
    {
    case SlipGeometry::kFCC:
      fSlipGeometry = new FCCGeometry(kSlipFCC);
      break;

    case SlipGeometry::kBCC:
      //fSlipGeometry = new BCCGeometry(kSlipBCC);
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: BCC not implemented");
      break;

    case SlipGeometry::kHCP:
      //fSlipGeometry = new HCPGeometry(kSlipHCP);
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: HCP not implemented");
      break;

    default:
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: Bad fCrystalType");
    }
  if (!fSlipGeometry) throwMemoryError("PolyCrystalMatT::SetSlipSystems");

  // set number of slip systems
  fNumSlip = fSlipGeometry->NumSlip();

  // allocate Schmidt tensor in crystal coords
  fZc.Allocate(fNumSlip);
  for (int i = 0; i < fNumSlip; i++)
    fZc[i].Allocate(3,3);

  // copy Schmidt tensor in crystal coords
  fZc = fSlipGeometry->GetSchmidtTensor();
}

void PolyCrystalMatT::SetLatticeOrientation()
{
  // input code to distribute ODF at each IP/ELEM
  fInput >> fODFCode;

  // incr to output texture
  fInput >> fODFOutInc;

  // read orientation data
  fLatticeOrient = new LatticeOrient(*this);
  if (!fLatticeOrient) throwMemoryError("PolyCrystalMatT::SetLatticeOrientation");

  // number of elements and integration points
  int numelem = NumElements();
  int numint = NumIP();

  // allocate array for euler angles at integration point
  fangles.Allocate(fNumGrain);
  for (int i = 0; i < fNumGrain; i++)
    fangles[i].Allocate(3);

  // allocate array to hold crystal orientations
  fEuler.Allocate(numelem);
  for (int i = 0; i< numelem; i++)
    {
      fEuler[i].Allocate(numint, fNumGrain);
      for (int j = 0; j < numint; j++)
	for (int k = 0; k < fNumGrain; k++)
	  fEuler[i](j,k).Allocate(3);
    }

  // assign orientation angles to each IP/ELEM
  fLatticeOrient->AssignEulerAngles(fODFCode, numelem, numint,
				    fNumGrain, fEuler);
}

void PolyCrystalMatT::SetCrystalElasticity()
{
  // input elasticity code
  fInput >> fElastCode;

  // select crystal elasticity type
  switch(fElastCode)
    {
    case CrystalElasticity::kIsoCrysElast:
      fElasticity = new IsotropicCrystalElast(*this);
      break;

    case CrystalElasticity::kCubCrysElast:
      fElasticity = new CubicCrystalElast(*this);
      break;

    default:
      throwRunTimeError("PolyCrystalMatT::SetCrystalElasticity: Bad fElastCode");
      break;
    }

  if (!fElasticity) throwMemoryError("PolyCrystalMatT::SetCrystalElasticity");

  // get elasticity constants (fmu, lambda, beta)
  fElasticity->ElasticityProps(fMatProp);
}

void PolyCrystalMatT::SetConstitutiveSolver()
{
  // input solver code
  fInput >> fSolverCode;

  // select constitutive solver
  switch(fSolverCode)
    {
      // Newton's method + line search
    case NLCSolver::kNLCSolver_LS: 
      fSolver = new NLCSolver_LS(fNumSlip);
      break;
          
      // Newton's method + trust region
    case NLCSolver::kNLCSolver_TR:
      //fSolver = new NLCSolver_TR(fNumSlip);
      throwRunTimeError("PolyCrystalMatT::SetConstitutiveSolver: TR not implemented");
      break; 

    default:
      throwRunTimeError("PolyCrystalMatT::SetConstitutiveSolver: Bad fSolverCode");
    }
  if (!fSolver) throwMemoryError("PolyCrystalMatT::SetConstitutiveSolver");

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
