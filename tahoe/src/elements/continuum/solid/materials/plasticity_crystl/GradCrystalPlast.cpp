/*
  File: GradCrystalPlast.cpp
*/

#include "GradCrystalPlast.h"
#include "SlipGeometry.h"
#include "LatticeOrient.h"
#include "NLCSolver.h"
#include "PowerLawIKinetics.h"
#include "PowerLawIIKinetics.h"
#include "VoceGradHardening.h"
#include "Utils.h"

#include "FEManagerT.h"
#include "ElementCardT.h"
#include "ifstreamT.h"
#include "ContinuumElementT.h"

//#include "ElasticT.h"
//DEV

//#include "ShapeFunctionT.h"
//DEV


/* spatial dimensions of the problem */
const int kNSD = 3;

const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

GradCrystalPlast::GradCrystalPlast(ifstreamT& in, const FiniteStrainT& element) :
  LocalCrystalPlast(in, element),  
  fLocInitX (ContinuumElement().InitialCoordinates()),
  fLocCurrX (LocalArrayT::kCurrCoords),
  fFeIP     (NumIP()),    
  fFeTrIP   (NumIP()),    
  fLDNa     (NumIP()),    
  fGDNa     (NumIP()),    
  fGradFe   (kNSD),
  fGradFeTr (kNSD),
  fGradU    (NumSD(),NumSD()), 
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

  // allocate space for ...
  // ... Fe values at integration points
  for (int i = 0; i < NumIP(); i++) {
    fFeIP[i].Allocate(kNSD,kNSD);
    fFeTrIP[i].Allocate(kNSD,kNSD);
  }

  // ... Fe values at nodal points
  fFeNodes.Allocate(fNumNodes);
  fFeTrNodes.Allocate(fNumNodes);
  for (int i = 0; i < fNumNodes; i++) {
    fFeNodes[i].Allocate(kNSD,kNSD);
    fFeTrNodes[i].Allocate(kNSD,kNSD);
  }

  // ... shape functions and their derivatives
  fLNa.Allocate(NumIP(), fNumNodes);
  for (int i = 0; i < NumIP(); i++) {
    fLDNa[i].Allocate(NumSD(), fNumNodes); 
    fGDNa[i].Allocate(NumSD(), fNumNodes); 
  }

  // ... nodal extrapolation matrix
  fNodalExtrap.Allocate(fNumNodes, NumIP());

  // ... spatial gradients of Fe (note: kNSD instead of NumSD())
  for (int i = 0; i < kNSD; i++) {
    fGradFe[i].Allocate(kNSD,kNSD);
    fGradFeTr[i].Allocate(kNSD,kNSD);
  }

  // ... current coordinates
  fLocCurrX.Allocate(fLocInitX.NumberOfNodes(), NumSD());

  // ... d(Curvature)/d(DGamma) (Moduli computation)
  fdKe.Allocate(NumIP());
  for (int i = 0; i < NumIP(); i++) {
    fdKe[i].Allocate(fNumSlip);
    for (int j = 0; j < fNumSlip; j++)
      fdKe[i][j].Allocate(kNSD,kNSD);
  }
}

GradCrystalPlast::~GradCrystalPlast() {} 

void GradCrystalPlast::Initialize()
{
  // inherited
  PolyCrystalMatT::Initialize();
  
  // shape funtions, their derivatives and extrap matrix in parent domain
  SetLocalShape(fLNa, fLDNa, fNodalExtrap);
}

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

  // compute crystal stresses
  if (fStatus == GlobalT::kFormRHS && CurrIP() == 0)
    {
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

	  // compute crystal Cauchy stress
	  CrystalS_ij();

	  // compute plastic deformation gradient
	  FPInverse();

	  // crystal moduli
	  fc_ijkl = 0.;
	  CrystalC_ijkl_Elastic();  // elastic part
	  CrystalC_ijkl_Plastic();  // plastic part
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
      
  // return crystal moduli
  return fc_ijkl;
}

void GradCrystalPlast::FormRHS(const dArrayT& dgamma, dArrayT& rhs)
{
  // incremental plastic deformation gradient ^ (-1)
  DeltaFPInverse(dgamma);

  // elastic right Cauchy-Green tensor 
  fCeBar.MultQTBQ(fDFpi, fCeBarTr);

  // resolved shear stress
  ResolveShearStress(fCeBar);

  // compute residual
  for (int i = 0; i < fNumSlip; i++)
    rhs[i] = dgamma[i] - fdt * fKinetics->Phi(fTau[i], i);
}

// form Jacobian 
void GradCrystalPlast::FormLHS(const dArrayT& dgamma, dMatrixT& lhs)
{
  // incremental plastic deformation gradient ^ (-1)
  //  DeltaFPInverse(dgamma);
  #pragma unused(dgamma)

  // elastic right Cauchy-Green tensor 
  //  fCeBar.MultQTBQ(fDFpi, fCeBarTr);

  // preliminary computations
  fdeltaI = 0.5 * fMatProp[1] * (fCeBar.Trace() - 3.) - fMatProp[0] + fMatProp[2];
  fCeBarTr.ToMatrix(fmatx1);
  fmatx4.MultATB(fDFpi, fmatx1); 

  // initialize jacobian
  lhs = 0.;

  // contribution of dTau/dDGamma to lhs(i,j)
  for (int i = 0; i < fNumSlip; i++)
    {
      farray[i].SetToScaled(-1., fZ[i]);

      // dTau/dCe
      dTaudCe(fZ[i], fP[i], fsymmatx1); 

      for (int j = 0; j < fNumSlip; j++)
	{
 	  // dCe/dDGamma
	  fmatx3.MultAB(fmatx4, fZ[j]);
	  fsymmatx2.Symmetrize(fmatx3);
	  fsymmatx2 *= -2.;

	  // dTau/dCe : dCe/dDGamma
          fsymmatx1.ToMatrix(fmatx1);
          fsymmatx2.ToMatrix(fmatx2);
          lhs(i,j) = dArrayT::Dot(fmatx1, fmatx2);
	}
    }

  // contribution of dPhi/dTau to lhs(i,i) 
  for (int i = 0; i < fNumSlip; i++)
    {
      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      for (int j = 0; j < fNumSlip; j++) lhs(i,j) *= -tmp;
      lhs(i,i) += 1.;
    }
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
  labels.Allocate(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void GradCrystalPlast::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
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
     cerr << " step # " << ContinuumElement().FEManager().StepNumber()
          << "    S_eq_avg = " 
          << sqrt(fsymmatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23 << endl; 

  // iteration counters
  output[1] = fIterCount;
  output[2] = fIterState;

  // compute euler angles
  const int& step = ContinuumElement().FEManager().StepNumber();
  const int& nsteps = ContinuumElement().FEManager().NumberOfSteps();

  if (fmod(step, fODFOutInc) == 0 || step == nsteps)
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
    fLatticeOrient->WriteTexture(elem, intpt, fNumGrain, step, fangles);
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

void GradCrystalPlast::InitializeCrystalVariables()
{
  // only one grain per integration point
  int igrn = 0;

  // initialize state at each element and ...
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);
      
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
  for (fIP = 0; fIP < NumIP(); fIP++) 
    {
      // deformation gradient at integration point IP
      //fFtot = DeformationGradient(fLocDisp); //DEV - deprecated
      Compute_Ftot_3D(fFtot, fIP);

      // fetch crystal data
      LoadCrystalData(element, fIP, igrn);
      
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
    }

  // iterate for the state in the element (all IP's at once)
  bool stateConverged = false;
  int iterState;
  for (iterState = 0; iterState < fMaxIterState; iterState++)
    {
      // compute dgamma and Fe at each integration point
      for (int intpt = 0; intpt < NumIP(); intpt++) 
	{
	  // fetch crystal data
	  LoadCrystalData(element, intpt, igrn);

	  // trial elastic right Cauchy-Green tensor
	  fCeBarTr.MultATA(fFeTr);

	  // Schmidt tensors and shear back stress
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fP[i].Symmetrize(fZ[i]);
	      fHardening->fTauKin[i] = dArrayT::Dot(fXe, fZ[i]);
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
	}
      
      // elastic curvature
      LatticeCurvature(element, igrn);

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

	  // solve for DDss and compute hardness
	  fHardening->ImplicitSolveHard();

	  // norm for hardness
	  fnormHard[intpt] = fHardening->Magnitude();
	}

      // check state convergence
      stateConverged = CheckConvergenceOfState();
      if (stateConverged) break;

      // reset norms
      fnormDGam0 = fnormDGam;
      fnormHard0 = fnormHard;
    }

  // check if did not converge in max iterations
  if (!stateConverged)
    throwRunTimeError("GradCrystalPlast::SolveCrystalState: Didn't converge in MaxIters");

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
    //    throw eGeneralFail;	
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

//DEV - deprecated
#if 0
const dMatrixT& GradCrystalPlast::DeformationGradient(const LocalArrayT& disp)
{
  ShapeFunction().GradU(disp, fGradU, fIP);
  return  FDContinuumT::F(fGradU);
}
#endif

/* PRIVATE MEMBER FUNCTIONS */

bool GradCrystalPlast::CheckConvergenceOfState() const
{
  int i = 0;
  bool test = true;
  while (i < NumIP() && test)
    {
      // test = ( fabs(fnormDGam[i]-fnormDGam0[i]) < fTolerState*fadots &&
      //          fabs(fnormHard[i]-fnormHard0[i]) < fTolerState*ftausi );
      test = ( fabs(fnormDGam[i]-fnormDGam0[i]) < fTolerState*fnormDGam0[i] &&
	       fabs(fnormHard[i]-fnormHard0[i]) < fTolerState*fnormHard0[i] );
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
  int blockToKe   = 4*kNSD*kNSD + fNumHard + fNumSlip;

  // current elastic curvature
  fKe.Set (kNSD,kNSD, &d_array[dex += blockToKe]);
}

void GradCrystalPlast::ShapeFunctionDeriv()
{
 // current coordinates of element nodes
  fLocCurrX.SetToCombination(1., fLocInitX, 1., fLocDisp);

  // derivarives dNa/dXcurr at integration points
  ComputeGDNa(fLocCurrX, fLDNa, fGDNa);
}

// dislocation tensor A2e in Bbar configuration
void GradCrystalPlast::LatticeCurvature(ElementCardT& element, int igrn)
{
  // extrapolate ipvalues of Fe to nodal points
  Extrapolate(fFeIP, fFeNodes, fNodalExtrap);

  for(int intpt = 0; intpt < NumIP(); intpt++)
    {
      // gradient of Fe
      ComputeGradFe(fGDNa[intpt], fFeNodes, fGradFe);

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
    }
}

void GradCrystalPlast::dKedDGamma(ElementCardT& element)
{
  // extrapolate ipvalues of Fe and FeTr to nodal points
  Extrapolate(fFeIP, fFeNodes, fNodalExtrap);
  Extrapolate(fFeTrIP, fFeTrNodes, fNodalExtrap);

  // compute derivative of curvature wrt DGamma
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // spatial gradients of Fe and FeTr at intpt
      ComputeGradFe(fGDNa[intpt], fFeNodes, fGradFe);
      ComputeGradFe(fGDNa[intpt], fFeTrNodes, fGradFeTr);

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


// dFe/dx at integration point
void GradCrystalPlast::ComputeGradFe(const dArray2DT& GDNa, const ArrayT<dMatrixT>& nodal,
                                     ArrayT<dMatrixT>& grad)
{
  // intialize spatial gradient d()/dx
  for (int i = 0; i < kNSD; i++) grad[i] = 0.;

  // compute dFe/dx at integration point
  if (NumSD() == 2)
    {
      double* dx1 = GDNa(0);
      double* dx2 = GDNa(1);
      
      for (int i = 0; i < fNumNodes; i++)
	{
	  grad[0].AddScaled(*dx1, nodal[i]);
	  grad[1].AddScaled(*dx2, nodal[i]);
	  
	  dx1++; dx2++;
	}
    }
  else // (NumSD() ==3)
    {
      double* dx1 = GDNa(0);
      double* dx2 = GDNa(1);
      double* dx3 = GDNa(2);
      
      for (int i = 0; i < fNumNodes; i++)
	{
	  grad[0].AddScaled(*dx1, nodal[i]);
	  grad[1].AddScaled(*dx2, nodal[i]);
	  grad[2].AddScaled(*dx3, nodal[i]);
	  
	  dx1++; dx2++; dx3++;
	}
    }
}

void GradCrystalPlast::CrystalC_ijkl_Plastic()
{
  // compute inverse of fDFpi
  fDFp.Inverse(fDFpi);

  // shear stress on slip systems
  ResolveShearStress(fCeBar);
  
  // SECOND TERM : (Fe_tr) (dDGamma/dCe_tr) (Fe_tr)^T
  FormLHS(fDGamma, fLHS);
  //AddToLHS(fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  for (int i = 0; i < fNumSlip; i++)
    {
      dTaudCe(fZ[i], fP[i], fsymmatx1);
      fA[i].MultQBQT(fFe, fsymmatx1);

      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      fA[i] *= -tmp;
    }

  for (int i = 0; i < fNumSlip; i++) {
    fB[i] = 0.;
    for (int j = 0; j < fNumSlip; j++)
      fB[i].AddScaled(-Jaci(i,j), fA[j]);
  }

  // FIRST TERM : (Fe_tr) (dS/dDGamma) (Fe_tr)^T
  for (int i = 0; i < fNumSlip; i++)
    {
      dMatrixT& matx = farray[i];      // dDFpi/dDGamma from FormLHS()
      fmatx4.MultAB(fDFp, matx);
      fmatx4 *= 2.;
      fsymmatx1.Symmetrize(fmatx4);
      dTaudCe(fmatx4, fsymmatx1, fsymmatx2);
      fA[i].MultQBQT(fFe, fsymmatx2);
    }
	
  // PLASTIC CONTRIBUTION TO MODULI : 2 * Sum(fA[i] (x) fB[i])
  for (int i = 0; i < fNumSlip; i++)
    {
      fRank4.Outer(fA[i], fB[i]);
      fc_ijkl.AddScaled(2., fRank4);
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


/*to be added to ShapeFunction class:

void ShapeFunctionT::GradU(const LocalArrayT& nodal, dMatrix& gradU, 
			   int IPnumber) const
{
  fDomain->Jacobian(nodal, (*pDNaU)[IPnumber], grad_U);
}
*/
