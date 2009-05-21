/* $Id: PolyCrystalMatT.cpp,v 1.19 2009-05-21 22:30:27 tdnguye Exp $ */
#include "PolyCrystalMatT.h"
#include "CrystalElasticity.h"
#include "SlipGeometry.h"
#include "LatticeOrient.h"
#include "NLCSolver.h"
#include "NLCSolver_LS.h"
#include "SlipKinetics.h"
#include "SlipHardening.h"
#include "Utils.h"
#include "ElementBlockDataT.h"

#include "FSMatSupportT.h"
#include "FiniteStrainT.h"
#include "StringT.h"

using namespace Tahoe;

/* number of elastic material properties : isotropic and cubic */
const int kNumMatProp = 3;

/* initialization flag value */
const int kIsInit = 1;

/* spatial dimensions of the problem */
const int kNSD = 3;

PolyCrystalMatT::PolyCrystalMatT(void) :
	ParameterInterfaceT("polycrystal_material"),
//  FDHookeanMatT(in, support),
  fdt           (0.0),
  fLocLastDisp  (NULL),
  fLocDisp      (NULL),
  fSlipGeometry (NULL),
  fLatticeOrient(NULL),
  fElasticity   (NULL),
  fKinetics     (NULL),
  fHardening    (NULL),
  fSolver       (NULL),
  fSolverPtr    (new SolverWrapperPoly(*this))
{

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

bool PolyCrystalMatT::NeedsPointInitialization() const { return true; }
void PolyCrystalMatT::PointInitialize(void)
{
	// only execute for first integration point
	if (CurrIP() == 0) {
	
		// determine storage size
		int i_size = NumIP();
		int d_size = NumVariablesPerElement();
	
		// current element
		ElementCardT& element = CurrentElement();
		
		// allocate space
		element.Dimension(i_size, d_size);
		element.IntegerData() = kIsInit;
		element.DoubleData()  = 0.0;

		// initialize state variables in all elements
		InitializeCrystalVariables(element);
	}
}

bool PolyCrystalMatT::NeedLastDisp() const { return true; }

/* describe the parameters needed by the interface */
void PolyCrystalMatT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FDHookeanMatT::DefineParameters(list);
	
	/* external parameters file */
	list.AddParameter(ParameterT::String, "parameter_file");
}

/* accept parameter list */
void PolyCrystalMatT::TakeParameterList(const ParameterListT& list)
{
	fLocLastDisp = &(FSMatSupport().FiniteStrain()->LastDisplacements());
	fLocDisp = &(FSMatSupport().FiniteStrain()->Displacements());

	/* dimension work space */
	fMatProp.Dimension(kNumMatProp);
	fFtot_n.Dimension(kNSD,kNSD);
	fFtot.Dimension(kNSD,kNSD);
	fFt.Dimension(kNSD,kNSD);
	fRotMat.Dimension(kNSD,kNSD);
	fs_ij.Dimension(kNSD);
	fsavg_ij.Dimension(kNSD);
	fc_ijkl.Dimension(dSymMatrixT::NumValues(kNSD));
	fcavg_ijkl.Dimension(dSymMatrixT::NumValues(kNSD));

	/* external file name */
	StringT filename = list.GetParameter("parameter_file");
	filename.Prepend(fstreamT::Root());
	fInput.set_marker('#');
	OpenExternal(fInput, filename, "PolyCrystalMatT data");

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

	// set slip system hardening law (set before slip kinetics!!!)
	SetSlipHardening();

	// set kinetics of slip
	SetSlipKinetics();
	
	// close file
	fInput.close();

	/* inherited - need to be called after SetCrystalElasticity */
	FDHookeanMatT::TakeParameterList(list);

#if 0
	// output stream
	ofstreamT& out = MaterialSupport().Output();

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
#endif
}

/* set (material) tangent modulus */
void PolyCrystalMatT::SetModulus(dMatrixT& modulus)
{
	fElasticity->ComputeModuli(modulus);
}

void PolyCrystalMatT::SetSlipSystems()
{
  // read crystal structure
  fInput >> fCrystalType;

  // construct slip system geometry
  switch(fCrystalType)
    {
    case SlipGeometry::kFCC:
      fSlipGeometry = new FCCGeometry();
      break;

    case SlipGeometry::kBCC:
      //fSlipGeometry = new BCCGeometry();
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: BCC not implemented");
      break;

    case SlipGeometry::kHCP:
      //fSlipGeometry = new HCPGeometry();
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: HCP not implemented");
      break;

    case SlipGeometry::kFCC24:
      fSlipGeometry = new FCCGeometry24();
      break;

    default:
      throwRunTimeError("PolyCrystalMatT::SetSlipSystems: Bad fCrystalType");
    }
  if (!fSlipGeometry) throwMemoryError("PolyCrystalMatT::SetSlipSystems");

  // set number of slip systems
  fNumSlip = fSlipGeometry->NumSlip();

  // allocate space for Schmidt tensor in crystal/sample coords
  fZc.Dimension(fNumSlip);
  fZ.Dimension(fNumSlip);
  for (int i = 0; i < fNumSlip; i++) {
    fZc[i].Dimension(3,3);
    fZ[i].Dimension(3,3);
  }

  // allocate space for slip vectors in crystal/sample coords
  fSlipSc.Dimension(fNumSlip);
  fSlipMc.Dimension(fNumSlip);
  fSlipS.Dimension(fNumSlip);
  fSlipM.Dimension(fNumSlip);
  for (int i = 0; i < fNumSlip; i++) {
    fSlipSc[i].Dimension(3);
    fSlipMc[i].Dimension(3);
    fSlipS[i].Dimension(3);
    fSlipM[i].Dimension(3);
  }

  // allocate space for slip shearing rate and resolve shear stress
  fDGamma.Dimension(fNumSlip);
  fGamma.Dimension(fNumSlip);
  fGamma = 0.0;
  fTau.Dimension(fNumSlip);

  // copy Schmidt tensor in crystal coords
  fZc = fSlipGeometry->GetSchmidtTensor();

  // copy slip vectors in crystal coords
  fSlipSc = fSlipGeometry->GetSlipVectorS();
  fSlipMc = fSlipGeometry->GetSlipVectorM();
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
  fangles.Dimension(fNumGrain);
  for (int i = 0; i < fNumGrain; i++)
    fangles[i].Dimension(3);

  // allocate array to hold crystal orientations
  fEuler.Dimension(numelem);
  for (int i = 0; i< numelem; i++)
    {
      fEuler[i].Dimension(numint, fNumGrain);
      for (int j = 0; j < numint; j++)
	for (int k = 0; k < fNumGrain; k++)
	  fEuler[i](j,k).Dimension(3);
    }

  
  // assign orientation angles to each elemgroup
  if (fODFCode == 5)
  {
	ArrayT<StringT> block_id;
	FSMatSupport().FiniteStrain()->ElementBlockIDs(block_id);
	int num_block = block_id.Length();
	for (int i = 0; i< num_block; i++)
	{
		const ElementBlockDataT& block_data = FSMatSupport().FiniteStrain()->BlockData(block_id[i]);
		int num_elem = block_data.Dimension();
		int start_elem = block_data.StartNumber();
		fLatticeOrient->AssignEulerAngles_block(num_elem, start_elem, numint, fNumGrain, fEuler);
	}
  }
  else
  // assign orientation angles to each IP/ELEM
	fLatticeOrient->AssignEulerAngles(fODFCode, numelem, numint, fNumGrain, fEuler);
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
  fdt = ContinuumElement().ElementSupport().TimeStep();
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
          if (totSubIncrs > 128) {
             cout << " ... in PolyCrystalMatT::SolveCrystalState: totSubIncrs > 128 \n";
             cout << " ...... will throw 'EBadJacobianDet' to force dtime decrease \n";
             throw ExceptionT::kBadJacobianDet;
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
      fdt = ContinuumElement().ElementSupport().TimeStep() * tmp;

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
		F_3D =  F_total();
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F_total());    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw ExceptionT::kGeneralFail;
}

void PolyCrystalMatT::Compute_Ftot_3D(dMatrixT& F_3D, int ip) const
{
	int nsd = NumSD();
	if (nsd == 3)
		F_3D =  F_total(ip);
	else if (nsd == 2)
	{
		// expand total deformation gradient: 2D -> 3D (plane strain)
		F_3D.Rank2ExpandFrom2D(F_total(ip));    // fFtot or fFtot_n
		F_3D(2, 2) = 1.0;
	}
	else 
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
		throw ExceptionT::kGeneralFail;
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
