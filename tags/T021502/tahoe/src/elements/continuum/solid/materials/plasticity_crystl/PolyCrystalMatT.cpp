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
  fFt           (kNSD,kNSD),
  fRotMat       (kNSD,kNSD),
  fs_ij         (kNSD),
  fsavg_ij      (kNSD),
  fc_ijkl       (dSymMatrixT::NumValues(kNSD)),
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

  // SetConstitutiveSolver moved to Initialize()
  // it allows to initialize NLCSolver with adequate # variables
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
  // set nonlinear constitutive solver 
  SetConstitutiveSolver();

  // read data for state iteration
  fInput >> fMaxIterState;
  fInput >> fTolerState;

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

  // allocate space for Schmidt tensor in crystal/sample coords
  fZc.Allocate(fNumSlip);
  fZ.Allocate(fNumSlip);
  for (int i = 0; i < fNumSlip; i++) {
    fZc[i].Allocate(3,3);
    fZ[i].Allocate(3,3);
  }

  // allocate space for slip shearing rate and resolve shear stress
  fDGamma.Allocate(fNumSlip);
  fTau.Allocate(fNumSlip);

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

  // number of variables in nonlinear solver
  int numvar = NumberOfUnknowns();

  // select constitutive solver
  switch(fSolverCode)
    {
      // Newton's method + line search
    case NLCSolver::kNLCSolver_LS: 
      fSolver = new NLCSolver_LS(numvar);
      break;
          
      // Newton's method + trust region
    case NLCSolver::kNLCSolver_TR:
      //fSolver = new NLCSolver_TR(numvar);
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

/* general driver to solve crystal state using subincrementation */
void PolyCrystalMatT::SolveCrystalState()
{
  // flag to track convergence of crystal state
  bool stateConverged;

  // counters for subincrementation 
  int subIncr = 1;
  int totSubIncrs = 1;

  // time step and deformation gradient
  fdt = ContinuumElement().FEManager().TimeStep();
  fFt = fFtot;

  // iterate to compute crystal state
  IterateOnCrystalState(stateConverged, subIncr);
 
  // if converged -> return; else -> do subincrementation
  if (stateConverged) return;
  
  // loop for subincrementation procedure
  for(;;)
    {
      // increase number of subincrements if not converged
      if (!stateConverged) 
	{
	  subIncr = 2 * subIncr - 1;
	  totSubIncrs = 2 * totSubIncrs;
          //if (totSubIncrs > pow(2,10)) throwRunTimeError
          //   ("PolyCrystalMatT::SolveCrystalState: totSubIncrs > 2^10");
          if (totSubIncrs > pow(2,10)) {
             cout << "PolyCrystalMatT::SolveCrystalState: totSubIncrs > 2^10 \n";
             cout << "  **will throw 'EBadJacobianDet' to force dtime decrease** \n";
             throw eBadJacobianDet;
          }
          if (subIncr > 1) RestoreSavedSolution();
	}

      // if converged, adjust subincrements
      else if (subIncr < totSubIncrs)
	{
	  if ((subIncr/2*2) == subIncr) 
	    {
	      subIncr = subIncr / 2 + 1;
	      totSubIncrs = totSubIncrs / 2;
	    }
	  else
	    subIncr = subIncr + 1;
	}

      // succesfull return for subincrementation
      else
	break;

      // time step
      double tmp = (float)subIncr / (float)totSubIncrs;
      fdt = ContinuumElement().FEManager().TimeStep() * tmp;

      // current deformation gradient
      fFt.SetToCombination( (1.-tmp), fFtot_n, tmp, fFtot );

      // save current converged solution before getting next solution 
      if (subIncr > 1 && stateConverged) {
         SaveCurrentSolution();
      }

      // iterate to compute crystal state
      IterateOnCrystalState(stateConverged, subIncr);
    }
}

/* compute 3D deformation gradient */
void PolyCrystalMatT::Compute_Ftot_3D(dMatrixT& F_3D) const
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

void PolyCrystalMatT::Compute_Ftot_3D(dMatrixT& F_3D, int ip) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D =  F(ip);
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F(ip));    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw eGeneralFail;
}

void PolyCrystalMatT::Compute_Ftot_last_3D(dMatrixT& F_3D) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D =  F_total_last();
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F_total_last());    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw eGeneralFail;
}

void PolyCrystalMatT::Compute_Ftot_last_3D(dMatrixT& F_3D, int ip) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D =  F_total_last(ip);
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F_total_last(ip));    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw eGeneralFail;
}

/* 4th order tensor transformation: Co_ijkl = F_iI F_jJ F_kK f_lL Ci_IJKL */
void PolyCrystalMatT::FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F)
{
  int nsd = F.Rows();
  dSymMatrixT coltemp;
  dMatrixT outer(nsd);
  dMatrixT transform(dSymMatrixT::NumValues(nsd));

  // compute tranformation matrix
  dArrayT Fi1(nsd,F(0));
  dArrayT Fi2(nsd,F(1));
  dArrayT Fi3(nsd,F(2));

  coltemp.Set(nsd,transform(0));
  outer.Outer(Fi1,Fi1);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(1));
  outer.Outer(Fi2,Fi2);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(2));
  outer.Outer(Fi3,Fi3);
  coltemp.FromMatrix(outer);

  coltemp.Set(nsd,transform(3));
  outer.Outer(Fi2,Fi3);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;

  coltemp.Set(nsd,transform(4));
  outer.Outer(Fi1,Fi3);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;

  coltemp.Set(nsd,transform(5));
  outer.Outer(Fi1,Fi2);
  outer.Symmetrize();
  coltemp.FromMatrix(outer);
  coltemp *= 2.0;
	
  // compute transformed tensor
  Co.MultQBQT(transform, Ci);
}
