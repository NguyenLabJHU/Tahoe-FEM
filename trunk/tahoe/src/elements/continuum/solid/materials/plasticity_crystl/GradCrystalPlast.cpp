/* $Id: GradCrystalPlast.cpp,v 1.11 2003-12-18 22:22:01 ebmarin Exp $ */
#include "GradCrystalPlast.h"
#include "SlipGeometry.h"
#include "LatticeOrient.h"
#include "NLCSolver.h"
#include "PowerLawIKinetics.h"
#include "PowerLawIIKinetics.h"
#include "VoceGradHardening.h"
#include "Utils.h"

#include "ElementCardT.h"
#include "ifstreamT.h"
#include "ContinuumElementT.h" //needed for ip coordinates

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 3;
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

/* to debug */
const bool XTAL_MESSAGES = false;
const int ELprnt = 0;

GradCrystalPlast::GradCrystalPlast(ifstreamT& in, const FSMatSupportT& support) :
  LocalCrystalPlast(in, support),  
  fLocInitX (ContinuumElement().InitialCoordinates()),
  fLocCurrX (LocalArrayT::kCurrCoords),
  fGradTool (NULL),
  fFeIP     (NumIP()),    
  fFeTrIP   (NumIP()),    
  fGradFe   (kNSD),
  fGradFeTr (kNSD),
  fKe_n     (kNSD,kNSD),  
  fKe       (kNSD,kNSD),  
  fXe       (kNSD,kNSD),  
  fnormDGam0(NumIP()),
  fnormHard0(NumIP()),
  fnormDGam (NumIP()),
  fnormHard (NumIP())
{
  // check number of grains
  if (fNumGrain != 1) 
    throwRunTimeError("GradCrystalPlast::GradCrystalPlast: NumGrain != 1");

  // number of element vertex nodes for gradient evaluation
  if (NumSD() == 2) 
    fNumNodes = 4;
  else if (NumSD() == 3) 
    fNumNodes = 8;
  else 
    throwRunTimeError("GradCrystalPlast::GradCrystalPlast: NumSD != 2 or 3");

  // create Gradient Tool object
  fGradTool = new GradientTools(NumIP(), fNumNodes, NumSD());

  // allocate space for ...
  // ... Fe values at integration points
  for (int i = 0; i < NumIP(); i++) {
    fFeIP[i].Dimension(kNSD,kNSD);
    fFeTrIP[i].Dimension(kNSD,kNSD);
  }

  // ... Fe values at nodal points
  fFeNodes.Dimension(fNumNodes);
  fFeTrNodes.Dimension(fNumNodes);
  for (int i = 0; i < fNumNodes; i++) {
    fFeNodes[i].Dimension(kNSD,kNSD);
    fFeTrNodes[i].Dimension(kNSD,kNSD);
  }

  // ... spatial gradients of Fe (note: kNSD instead of NumSD())
  for (int i = 0; i < kNSD; i++) {
    fGradFe[i].Dimension(kNSD,kNSD);
    fGradFeTr[i].Dimension(kNSD,kNSD);
  }

  // ... current coordinates
  fLocCurrX.Dimension(fLocInitX.NumberOfNodes(), NumSD());

  // ... d(Curvature)/d(DGamma) (Moduli computation)
  fdKe.Dimension(NumIP());
  for (int i = 0; i < NumIP(); i++) {
    fdKe[i].Dimension(fNumSlip);
    for (int j = 0; j < fNumSlip; j++)
      fdKe[i][j].Dimension(kNSD,kNSD);
  }
}

GradCrystalPlast::~GradCrystalPlast() {} 

int GradCrystalPlast::NumVariablesPerElement()
{
  int d_size = 0;
  int strdim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                        // fRotMat (const)
  d_size += kNSD*kNSD;                        // fFpi_n   (t_n)
  d_size += kNSD*kNSD;                        // fKe_n    (t_n)
  d_size += fNumSlip;                         // fDGamma_n(t_n) 
  d_size += kNSD*kNSD;                        // fFpi     (t)
  d_size += kNSD*kNSD;                        // fKe      (t)
  d_size += fNumSlip;                         // fDGamma  (t)
  d_size += kNSD*kNSD;                        // fFeTr    (t)
  d_size += kNSD*kNSD;                        // fXe      (t)
  d_size += strdim;                           // fs_ij    (t)
  d_size += strdim*strdim;                    // fc_ijkl  (t)
  d_size += fHardening->NumberOfVariables();  // hard vars at t and t_n

  // total # variables per element
  d_size *= fNumGrain * NumIP();

  return d_size;
}

const dSymMatrixT& GradCrystalPlast::s_ij()
{
  // fetch current element
  ElementCardT & element = CurrentElement();

  // only one grain per integration point
  int igrn = 0;

  // time step
  fdt = fFSMatSupport.TimeStep();

  // compute crystal stresses
  if (fFSMatSupport.RunState() == GlobalT::kFormRHS && CurrIP() == 0)
    {
      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt)
         cout << " elem # " << CurrElementNumber() << endl;

      // state at each integration point (all at once)
      SolveCrystalState();

      // crystal stresses and moduli (all at once)
      for (fIP = 0; fIP < NumIP(); fIP++)
	{
	  // fetch crystal data
	  LoadCrystalData(element, fIP, igrn);
	  
	  // Schmidt tensor and density of GND's
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fHardening->fTauKin[i] = dArrayT::Dot(fXe, fZ[i]);
	    }

          // elasticity matrix in Bbar configuration
          if (!fElasticity->IsIsotropic())
             {
                fElasticity->ComputeModuli(fRank4);
                FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
             }
          else
                fElasticity->ComputeModuli(fcBar_ijkl);

	  // compute crystal Cauchy stress
	  CrystalS_ij();

	  // compute plastic deformation gradient
	  FPInverse();

	  // crystal moduli
	  fc_ijkl = 0.;
	  CrystalC_ijkl();

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt)
             cout << " IP # " << fIP << endl
                  << "s_ij \n  " << fs_ij << endl
                  << "c_ijkl \n" << fc_ijkl << endl;
	}
    }

  // fetch crystal data
  LoadCrystalData(element, CurrIP(), igrn);

  // return crystal stress
  return fs_ij;
}

const dMatrixT& GradCrystalPlast::c_ijkl()
{
  // fetch current element and intpt
  ElementCardT & element = CurrentElement();
  int intpt = CurrIP();

  // only one grain per integration point
  int igrn = 0;

  // recover local data
  LoadCrystalData(element, intpt, igrn);
      
  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
       cout << " GradCrystalPlast::c_ijkl " << endl;
       cout << " IP # " << intpt << endl
            << "c_ijkl \n" << fc_ijkl << endl;
  }

  // return crystal moduli
  return fc_ijkl;
}

void GradCrystalPlast::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // only one grain per integration point
  int igrn = 0;

  // update state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadCrystalData(element, intpt, igrn);
      
      // update state
      fFpi_n    = fFpi;
      fKe_n     = fKe;
      fDGamma_n = fDGamma;
      fHardening->UpdateHistory();
    }
}

void GradCrystalPlast::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // only one grain per integration point
  int igrn = 0;

  // reset state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadCrystalData(element, intpt, igrn);
      
      // reset state
      fFpi    = fFpi_n;
      fKe     = fKe_n;
      fDGamma = fDGamma_n;
      fHardening->ResetHistory();
    }
}

int GradCrystalPlast::NumOutputVariables() const {return kNumOutput;}

void GradCrystalPlast::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void GradCrystalPlast::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int group = ContinuumElement().ElementGroupNumber();
  int elem  = CurrElementNumber();
  int intpt = CurrIP();

  // only one grain per integration point
  int igrn = 0;

  // fetch crystal data 
  LoadCrystalData(element, intpt, igrn);
      
  // Von Mises stress
  output[0] = sqrt(fsymmatx1.Deviatoric(fs_ij).ScalarProduct())/sqrt23;
  // cerr << " S_eq = " << output[0] << endl;

  // compute averaged equivalent stress
  if (elem == 0 && intpt == 0) fAvgStress = 0.0;
  fAvgStress.AddScaled(1./(NumIP()*NumElements()), fs_ij);
  if (elem == (NumElements()-1) && intpt == (NumIP()-1))
     cerr << " step # " << fFSMatSupport.StepNumber()
          << "    S_eq_avg = " 
          << sqrt(fsymmatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23 << endl; 

  // iteration counters
  output[1] = fIterCount;
  output[2] = fIterState;

  // compute euler angles
  int step   = fFSMatSupport.StepNumber();
  int nsteps = fFSMatSupport.NumberOfSteps();

  if (fmod(double(step), fODFOutInc) == 0 || step == nsteps)
  {
    // elastic deformation gradient
    DeltaFPInverse(fDGamma);
    fFe.MultAB(fFeTr, fDFpi);

    // texture: rotation tensor from fFe
    PolarDecomp();
    
    // texture: compute new crystal orientation matrix
    fmatx1.MultAB(fRe, fRotMat);

    // texture: compute Euler angles (in radians)
    dArrayT& angles = fangles[igrn];
    fLatticeOrient->RotMatrixToAngles(fmatx1, angles);

    // write euler angles at IP/ELE
    fLatticeOrient->WriteTexture(group, elem, intpt, fNumGrain, step, fangles);
  }
}

void GradCrystalPlast::Print(ostream& out) const
{
  // inherited
  PolyCrystalMatT::Print(out);

  // print slip kinetics data
  out << "    Crystal slip kinetics (gradient crystal plast.)\n";
  out << "       Kinetics law. . . . . . . . . . . . . . . = " << fKinEqnCode << "\n";
  fKinetics->Print(out);

  // print slip hardening data
  out << "    Crystal slip hardening (gradient crystal plast.)\n";
  out << "       Hardening law . . . . . . . . . . . . . . = " << fHardCode << "\n";
  fHardening->Print(out);
}

void GradCrystalPlast::PrintName(ostream& out) const
{
  // inherited
  PolyCrystalMatT::PrintName(out);

  // output model name
  out << "    gradient crystal plasticity equations\n";
  fSlipGeometry->PrintName(out);
  fKinetics->PrintName(out);
  fHardening->PrintName(out);
}

/* PROTECTED MEMBER FUNCTIONS */

void GradCrystalPlast::SetSlipKinetics()
{
  // read slip kinetics code model
  fInput >> fKinEqnCode;

  // select slip hardening law
  switch(fKinEqnCode)
    {
    case SlipKinetics::kPowLawI:         // standard power law (iso)
      fKinetics = new PowerLawIKinetics(*this);
      break;
      
    case SlipKinetics::kPowLawII:        // standard power law (iso + kin)
      fKinetics = new PowerLawIIKinetics(*this);
      break;

    default:
      throwRunTimeError("GradCrystalPlast::SetSlipKinetics: Bad fKinEqnCode");
    }
  if (!fKinetics) throwMemoryError("GradCrystalPlast::SetSlipKinetics");
}

void GradCrystalPlast::SetSlipHardening()
{
  // read hardening code model
  fInput >> fHardCode;

  // select slip hardening law
  switch(fHardCode)
    {
    case SlipHardening::kHard_G1:   // ss-gn dislocation type model
      fHardening = new VoceGradHardening(*this);
      break;
      
    default:
      throwRunTimeError("GradCrystalPlast::SetSlipHardening: Bad fHardCode");
    }
  if (!fHardening) throwMemoryError("GradCrystalPlast::SetSlipHardening");
}

void GradCrystalPlast::InitializeCrystalVariables(ElementCardT& element)
{
  // only one grain per integration point
  int igrn = 0;
  int elem = CurrElementNumber();

	// ... at each integration point
	for (int intpt = 0; intpt < NumIP(); intpt++)
	{
	  // fetch crystal data 
	  LoadCrystalData(element, intpt, igrn);
	  
	  // fetch euler angles
	  dArrayT& angles = fEuler[elem](intpt, igrn);
	  
	  // storage rotation matrix from Euler angles
	  fLatticeOrient->AnglesToRotMatrix(angles, fRotMat);
	  
	  // plastic deformation gradients
	  fFpi_n.Identity();
	  fFpi.Identity();

	  // elastic curvatures
	  fKe_n = 0.;
	  fKe   = 0.;
	  
	  // shear rates on slip systems
	  fDGamma_n = 0.;
	  fDGamma   = 0.;
	      
	  // elastic deformation gradient
	  fFeTr.Identity();
	      
	  // stress associated with lattice curvature
	  fXe = 0.;
	      
	  // crystal Cauchy stress and moduli
	  fs_ij = 0.;
	  fc_ijkl = 0.;

	  // hardening variables
	  fHardening->InitializeHardVariables();
	}
}

void GradCrystalPlast::LoadCrystalData(ElementCardT& element,
				       int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int strdim      = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 7*kNSD*kNSD + strdim + strdim*strdim + 2*fNumSlip
                   + fHardening->NumberOfVariables();
  int dex         = intpt*fNumGrain*blockPerGrn + igrain*blockPerGrn;

  fRotMat.Set    (kNSD,kNSD,     &d_array[dex             ]);   
  fFpi_n.Set     (kNSD,kNSD,     &d_array[dex += kNSD*kNSD]);     
  fKe_n.Set      (kNSD,kNSD,     &d_array[dex += kNSD*kNSD]);   
  fDGamma_n.Set  (fNumSlip,      &d_array[dex += kNSD*kNSD]);
  fFpi.Set       (kNSD,kNSD,     &d_array[dex += fNumSlip ]);     
  fKe.Set        (kNSD,kNSD,     &d_array[dex += kNSD*kNSD]);   
  fDGamma.Set    (fNumSlip,      &d_array[dex += kNSD*kNSD]);
  fFeTr.Set      (kNSD,kNSD,     &d_array[dex += fNumSlip ]);
  fXe.Set        (kNSD,kNSD,     &d_array[dex += kNSD*kNSD]);   
  fs_ij.Set      (kNSD,          &d_array[dex += kNSD*kNSD]);
  fc_ijkl.Set    (strdim,strdim, &d_array[dex += strdim   ]);
  fHardening->LoadHardData(strdim*strdim, dex, d_array);
}

void GradCrystalPlast::SolveCrystalState()
{
  // fetch current element
  ElementCardT& element = CurrentElement();

  // reset iteration counter to check NLCSolver
  fIterState = 0;
  fIterCount = 0;

  // only one grain per integration point
  int igrn = 0;

  // shape function derivatives to compute GradFe
  ShapeFunctionDeriv();

  // get reference to hardening material properties
  const dArrayT& prop = fHardening->MaterialProperties();

  // elastic predictor Fe and initial guess for (dgamma, hardness) 
  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt)
         cout << "\n\n***Elastic Predictor Phase" << endl;

  for (fIP = 0; fIP < NumIP(); fIP++) 
    {
      // deformation gradient at integration point IP
      //fFtot = DeformationGradient(fLocDisp); //DEV - deprecated
      Compute_Ftot_3D(fFtot, fIP);
      Compute_Ftot_last_3D(fFtot_n, fIP);

      // fetch crystal data
      LoadCrystalData(element, fIP, igrn);
      
      // elasticity matrix in Bbar configuration
      if (!fElasticity->IsIsotropic())
         {
            fElasticity->ComputeModuli(fRank4);
            FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
         }
      else
            fElasticity->ComputeModuli(fcBar_ijkl);

      // trial elastic deformation gradient
      fFeTr.MultAB(fFtot, fFpi_n);

      // bookeeping to compute Grad(FeTr) (for computing c_ijkl) 
      fFeTrIP[fIP] = fFeTr;

      // unsymmetric back stress fXe
      fXe.SetToScaled(prop[4]*prop[5]*fMatProp[0], fKe_n); 
      
      // Schmidt tensor and density of GND's
      for (int i = 0; i < fNumSlip; i++) 
	{
	  fZ[i].MultQBQT(fRotMat, fZc[i]);
	  fHardening->fTauKin[i] = dArrayT::Dot(fXe, fZ[i]);
	}

      // initial estimate for DGamma
      ForwardGradientEstimate();

      // initial guess for hardness
      fHardening->ExplicitUpdateHard();

      // initial norms for dgamma and hardness
      fnormDGam0[fIP] = fDGamma.Magnitude();
      fnormHard0[fIP] = fHardening->Magnitude();

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << " IP # " << fIP << "\n"
              << "    fFtot: \n" << fFtot << endl;
         cout << "    fFeTr: \n" << fFeTr << endl;
         cout << "    fKe_n: \n" << fKe_n << endl;
         cout << "    fXe: \n" << fXe << endl;
         cout << "    fTauKin : \n" << fHardening->fTauKin << endl;
         cout << "    fDGamma (Fwd): \n" << fDGamma << endl;
         cout << "    normDGam0 : " << fnormDGam0[fIP] << endl;
         cout << "    normHard0 : " << fnormHard0[fIP] << endl;
      }
    }

  // iterate for the state in the element (all IP's at once)
  bool stateConverged = false;
  int iterState;
  for (iterState = 0; iterState < fMaxIterState; iterState++)
    {
      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << "\n\n***Iterate for state" << endl;
         cout << "   iterState # " << iterState << endl;
      }

      // compute dgamma and Fe at each integration point
      for (int intpt = 0; intpt < NumIP(); intpt++) 
	{
	  // fetch crystal data
	  LoadCrystalData(element, intpt, igrn);

          // elasticity matrix in Bbar configuration
          if (!fElasticity->IsIsotropic())
             {
                fElasticity->ComputeModuli(fRank4);
                FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
             }
          else
                fElasticity->ComputeModuli(fcBar_ijkl);

	  // trial elastic right Cauchy-Green tensor
	  fCeBarTr.MultATA(fFeTr);

	  // Schmidt tensors and shear back stress
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fP[i].Symmetrize(fZ[i]);
	      fHardening->fTauKin[i] = dArrayT::Dot(fXe, fZ[i]);
	    }
	  
          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             cout << " IP # " << intpt << "\n";
             cout << "    fFeTr: \n" << fFeTr << endl;
             cout << "    fTauKin : \n" << fHardening->fTauKin << endl;
          }

	  // solve for shear deformation on slip systems
	  SolveForDGamma();
	  
	  // incremental plastic deformation gradient ^ (-1)
	  DeltaFPInverse(fDGamma);
	  
	  // Elastic Deformation Gradient
	  fFe.MultAB(fFeTr, fDFpi);

	  // bookeeping to compute lattice curvature
	  fFeIP[intpt] = fFe;

	  // norm for dgamma
	  fnormDGam[intpt] = fDGamma.Magnitude();

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             cout << "    fDGamma (after solve): \n" << fDGamma << endl;
             cout << "    fDFpi: \n" << fDFpi << endl;
             cout << "    fFe: \n" << fFe << endl;
             cout << "    normDGam : " << fnormDGam[intpt] << endl;
          }
	}
      
      // elastic curvature
      LatticeCurvature(element, igrn);

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) 
         cout << "\n\n***Compute Hardness" << endl;

      // compute hardness
      for (int intpt = 0; intpt < NumIP(); intpt++) 
	{
	  // fetch crystal data
	  LoadCrystalData(element, intpt, igrn);
	  
          // unsymmetric back stress fXe
          fXe.SetToScaled(prop[4]*prop[5]*fMatProp[0], fKe); 
      
	  // Schmidt tensor and density of GND's
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fHardening->fTauKin[i] = dArrayT::Dot(fXe, fZ[i]);
	    }

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             cout << " IP # " << intpt << "\n";
             cout << "    fXe: \n" << fXe << endl;
             cout << "    fTauKin : \n" << fHardening->fTauKin << endl;
          }
	  // solve for DDss and compute hardness
	  //fHardening->ImplicitSolveHard();
	  fHardening->ImplicitUpdateHard();

	  // norm for hardness
	  fnormHard[intpt] = fHardening->Magnitude();

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             //cout << " IP # " << intpt << "\n";
             //cout << "    fXe: \n" << fXe << endl;
             //cout << "    fTauKin : \n" << fHardening->fTauKin << endl;
             cout << "    normHard : " << fnormHard[intpt] << endl;
          }
	}

      // check state convergence
      stateConverged = CheckConvergenceOfState();

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << "\n\n***Convergence Check: \n"; 
         cout << " stateConverged: " << stateConverged << endl;
         for (int i=0; i<NumIP(); i++)
           cout << " IP # " << i  << endl
                << "     normHard0 : " << fnormHard0[i] << endl
                << "     normHard : "  << fnormHard[i] << endl
                << "     normDGam0 : " << fnormDGam0[i] << endl
                << "     normDGam : "  << fnormDGam[i] << endl
                << "     | diff normGam  | = " << fabs(fnormDGam[i]-fnormDGam0[i]) << endl
                << "     | diff normHard | = " << fabs(fnormHard[i]-fnormHard0[i]) << endl;
      }

      if (stateConverged) break;

      // reset norms
      fnormDGam0 = fnormDGam;
      fnormHard0 = fnormHard;
    }


  // check if did not converge in max iterations
  if (!stateConverged) {
    cout << " elem # " << CurrElementNumber() << endl;
    throwRunTimeError("GradCrystalPlast::SolveCrystalState: Didn't converge in MaxIters");
  }

  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt)
      cout << " \n\n STATE CONVERGED " << endl;

  // update iteration count for state
  fIterState = max(fIterState, iterState);

  // compute dKe/dDGamma (for Jacobian expression)
  dKedDGamma(element);
}

void GradCrystalPlast::SolveForDGamma()
{
  int ierr = 0;
  dArrayT dgamma_n = fDGamma;

  // solve for incremental shear strain
  fSolver->Solve(fSolverPtr, fDGamma, ierr);

  if (ierr == 1) {
    writeMessage("GradCrystalPlast::SolveForDGamma: Convergence problems");
    writeWarning("GradcrystalPlast::SolveForDGamma: Will set:  dgamma = dgamma_n");
    fDGamma = dgamma_n;
    //    throw ExceptionT::kGeneralFail;	
  }

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());
}

void GradCrystalPlast::DeltaFPInverse(const dArrayT& dgamma)
{
  // incremental plastic deformation gradient ^ (-1)
  fDFpi.Identity();
  for (int i = 0; i < fNumSlip; i++)
    {
      fDFpi.AddScaled(-dgamma[i], fZ[i]);
    }
}

/* PRIVATE MEMBER FUNCTIONS */

bool GradCrystalPlast::CheckConvergenceOfState() const
{
  int i = 0;
  bool test = true;
  while (i < NumIP() && test)
    {
      // test = ( fabs(fnormDGam[i]-fnormDGam0[i]) < fTolerState*fadots &&
      //          fabs(fnormHard[i]-fnormHard0[i]) < fTolerState*ftausi );
      test = ( fabs(fnormDGam[i]-fnormDGam0[i]) <= fTolerState*fnormDGam0[i] &&
	       fabs(fnormHard[i]-fnormHard0[i]) <= fTolerState*fnormHard0[i] );
      i++;
    }
  return test;
}

void GradCrystalPlast::LoadCrystalCurvature(ElementCardT& element, 
					    int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int strdim      = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 7*kNSD*kNSD + strdim + strdim*strdim + 2*fNumSlip
                   + fHardening->NumberOfVariables();
  int dex         = intpt*fNumGrain*blockPerGrn + igrain*blockPerGrn;
  int blockToKe   = 4*kNSD*kNSD + fNumSlip;

  // current elastic curvature
  fKe.Set (kNSD,kNSD, &d_array[dex += blockToKe]);
}

void GradCrystalPlast::ShapeFunctionDeriv()
{
 // current coordinates of element nodes
  fLocCurrX.SetToCombination(1., fLocInitX, 1., fLocDisp);

  // derivarives dNa/dXcurr at integration points
  fGradTool->ComputeGDNa(fLocCurrX);

  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
     cout << "***\n\nShapeFunctionDerivatives" << endl;
     cout << "\n fLocDisp  (#nodes,#sd) :\n " << fLocDisp  << endl;
     cout << "\n fLocCurrX (#nodes,#sd) :\n " << fLocCurrX << endl;
  //   for (int i=0; i<NumIP(); i++) 
  //      cout << " fGDNa for ip # " << i << endl << fGDNa[i] << endl;
  }
}

// dislocation tensor A2e in Bbar configuration
void GradCrystalPlast::LatticeCurvature(ElementCardT& element, int igrn)
{
  // extrapolate ipvalues of Fe to nodal points
  fGradTool->ExtrapolateTensor(fFeIP, fFeNodes);

  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
     cout << "\n\n***Lattice Curvature " << endl;
     for (int i=0; i<fNumNodes; i++) 
        cout << " FeNodes at node # " << i << endl << fFeNodes[i] << endl;   
  }

  for(int intpt = 0; intpt < NumIP(); intpt++)
    {
      // gradient of Fe
      fGradTool->GradientOfTensorAtIP(fFeNodes, fGradFe, intpt);

      // fetch crystal curvature
      LoadCrystalCurvature(element, intpt, igrn);

      // auxiliary matrices
      dMatrixT& Fe  = fFeIP[intpt];
      dMatrixT& Fei = fmatx1.Inverse(fFeIP[intpt]);

      // curvature tensor fKe (dislocation tensor ABar2e)
      fKe = 0;
      for (int I = 0; I < kNSD; I++)
	for (int j = 0; j < kNSD; j++)
	  {
	    dMatrixT& Fe_x = fGradFe[j];
	    for (int i = 0; i < kNSD; i++)
	      {
		fKe(0,I) += (-Fe_x(i,1)*Fe(j,2)+Fe_x(i,2)*Fe(j,1))*Fei(I,i);
		fKe(1,I) += ( Fe_x(i,0)*Fe(j,2)-Fe_x(i,2)*Fe(j,0))*Fei(I,i);
		fKe(2,I) += (-Fe_x(i,0)*Fe(j,1)+Fe_x(i,1)*Fe(j,0))*Fei(I,i);
	      }
	  }

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << " IP # " << intpt << "\n";
         for(int i=0; i<kNSD; i++)
            cout << "   fGradFe for dir " << i << endl << fGradFe[i] << endl;
         cout << "    Curvature tensor fKe: \n" << fKe << endl;
      }
    }
}

void GradCrystalPlast::dKedDGamma(ElementCardT& element)
{
  // extrapolate ipvalues of Fe and FeTr to nodal points
  fGradTool->ExtrapolateTensor(fFeIP, fFeNodes);
  fGradTool->ExtrapolateTensor(fFeTrIP, fFeTrNodes);

  // compute derivative of curvature wrt DGamma
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // spatial gradients of Fe and FeTr at intpt
      fGradTool->GradientOfTensorAtIP(fFeNodes, fGradFe, intpt);
      fGradTool->GradientOfTensorAtIP(fFeTrNodes, fGradFeTr, intpt);

      // auxiliar references
      dMatrixT& Fe    = fFeIP[intpt];
      dMatrixT& Fei   = fmatx1.Inverse(fFeIP[intpt]);
      dMatrixT& FeTr  = fFeTrIP[intpt];
      dMatrixT& FeiTr = fmatx2.Inverse(fFeTrIP[intpt]);

      // fetch crystal data
      LoadCrystalData(element, intpt, 0);    // igrn = 0

      for (int islip = 0; islip < fNumSlip; islip++)
	{
	  // Schmidt tensor
	  fmatx3.MultQBQT(fRotMat, fZc[islip]);
	  dMatrixT& Z = fmatx3;

	  // reference to dKedDGamma array
	  dMatrixT& dKe = fdKe[intpt][islip];
          dKe = 0.0;

	  // compute dKe/dDGamma at intpt
	  for (int I = 0; I < kNSD; I++)
	    for (int i = 0; i < kNSD; i++)
	      for (int l = 0; l < kNSD; l++)
		{
		  dMatrixT& Fe_x   = fGradFe[l];
		  dMatrixT& FeTr_x = fGradFeTr[l];
		  
		  for (int M = 0; M < kNSD; M++)
		    {
		      dKe(0,M) += (-(-Z(I,1))*Fe(l,2) + (-Z(I,2))*Fe(l,1)) * FeTr_x(i,I) * Fei(M,i)
			+ (-Fe_x(i,1)*(-Z(I,2)) + Fe_x(i,2)*(-Z(I,1))) * FeTr(l,I) * Fei(M,i)
			+ ( Fe_x(i,1)*Fe(l,2) - Fe_x(i,2)*Fe(l,1)) * (-1.) * Z(M,I) * FeiTr(I,i);
		      
		      dKe(1,M) += ((-Z(I,0))*Fe(l,2) - (-Z(I,2))*Fe(l,0)) * FeTr_x(i,I) * Fei(M,i)
			+ ( Fe_x(i,0)*(-Z(I,2)) - Fe_x(i,2)*(-Z(I,0))) * FeTr(l,I) * Fei(M,i)
			+ (-Fe_x(i,0)*Fe(l,2) + Fe_x(i,2)*Fe(l,0)) * (-1.) * Z(M,I) * FeiTr(I,i);
		      
		      dKe(2,M) += (-(-Z(I,0))*Fe(l,1) + (-Z(I,1))*Fe(l,0)) * FeTr_x(i,I) * Fei(M,i)
			+ (-Fe_x(i,0)*(-Z(I,1)) + Fe_x(i,1)*(-Z(I,0))) * FeTr(l,I) * Fei(M,i)
			+ ( Fe_x(i,0)*Fe(l,1) - Fe_x(i,1)*Fe(l,0)) * (-1.) * Z(M,I) * FeiTr(I,i);
		    }
		}
	}
    }
}

void GradCrystalPlast::AddToLHS(dMatrixT& lhs)
{
  // compute material constant
  const dArrayT& prop = fHardening->MaterialProperties();
  double cxbmu = prop[4] * prop[5] * fMatProp[0];

  // contribution of d(Phi)/d(DGamma)
  dArrayT& tauX = fHardening->fTauKin;
  for (int is = 0; is < fNumSlip; is++)
    {
      // compute derivative terms
      double dphidtaus = fdt * fKinetics->DPhiDIso(fTau[is], is);
      double dphidtaux = fdt * fKinetics->DPhiDKin(fTau[is], is);
      double dHddgamma = HardFuncDerivative(fDGamma[is],fHardening->fTauIso[is],tauX[is],1);
      double dHdtaus   = HardFuncDerivative(fDGamma[is],fHardening->fTauIso[is],tauX[is],2);
      double dHdtaux   = HardFuncDerivative(fDGamma[is],fHardening->fTauIso[is],tauX[is],3);
      double h = 1. - fdt * dHdtaus;

      for (int js = 0; js < fNumSlip; js++)
	{
	  // derivative d(Ke)/d(DGamma)
	  dMatrixT& dKe = fdKe[fIP][js];

	  // add -(d(Phi)/d(TauX)+dt/h*d(Phi)/d(TauS)*d(H)/dTauX) * d(TauX)/d(DGamma),
	  // where TauX = (cx*b*mu)*Ke:Z
	  lhs(is, js) -= (dphidtaux + fdt/h*dphidtaus*dHdtaux)*cxbmu*dArrayT::Dot(dKe, fZ[is]);

	  // add diagonal term -dt/h * d(Phi)/TauS * dH/dDGamma
	  if (is == js) lhs(is, js) -= fdt/h*dphidtaus*dHddgamma;
	}
    }
}

double GradCrystalPlast::HardFuncDerivative(double& dgamma, double& taus,
					       double& taux, int kcode)
{
  double dHdv;

  // fetch hardening material constants
  const dArrayT& prop = fHardening->MaterialProperties();

  // some constants
  double coeff1 = prop[3]*prop[0]/prop[4];
  double coeff2 = prop[0]/(prop[2]-prop[1]);

  // accumulated quantities
  double wkrate = 0.;
  double shrate = 0.;
  for (int i = 0; i < fNumSlip; i++) { 
    shrate += fabs(fDGamma[i]);
    wkrate += fabs(fHardening->fTauKin[i])*fabs(fDGamma[i]);
  }

  switch (kcode)
    {
    case 1:          // dTheta/dDGamma
      dHdv = ( coeff1/(taus-prop[1])*fabs(taux)
		+ coeff2*(prop[2]-taus) ) * fabs(dgamma)/dgamma;
      break;
      
    case 2:         // dTheta/dTauS
      dHdv = -coeff1/((taus-prop[1])*(taus-prop[1]))*wkrate - coeff2*shrate;
      break;

    case 3:         // dTheta/dTauX
      dHdv = coeff1/(taus-prop[1]) * fabs(dgamma) * fabs(taux)/taux;
      break;

    default:
      throwRunTimeError("GradCrystalPlast::HardFuncDerivative: Bad kcode");
    }

  return dHdv;
}
