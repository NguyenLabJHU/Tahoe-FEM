/*
  File: LocalCrystalPlast.cpp
*/

#include "LocalCrystalPlast.h"
#include "SlipGeometry.h"
#include "LatticeOrient.h"
#include "CrystalElasticity.h"
#include "NLCSolver.h"
#include "PowerLawIKinetics.h"
#include "PowerLawIIKinetics.h"
#include "HaasenKinetics.h"
#include "VoceHardening.h"
#include "HaasenHardening.h"
#include "ElementCardT.h"
#include "ifstreamT.h"
#include "Utils.h"

#include "ElasticT.h"
#include "FEManagerT.h"

/* spatial dimensions of the problem */
const int kNSD = 3;

const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

LocalCrystalPlast::LocalCrystalPlast(ifstreamT& in, const ElasticT& element) :
  PolyCrystalMatT(in, element),  

  // deformation gradient (continuation method)
  fFt   (kNSD,kNSD),

  // elastic deformation gradients
  fFeTr (kNSD,kNSD),
  fFe   (kNSD,kNSD),

  // plastic deformation gradients 
  fFpi_n (kNSD,kNSD),
  fFpi   (kNSD,kNSD),
  fDFp   (kNSD,kNSD),
  fDFpi  (kNSD,kNSD),

  // symmetric tensors in interm config	
  fCeBarTr (kNSD),
  fCeBar   (kNSD),
  fEeBar   (kNSD),
  fSBar    (kNSD),

  // crystal orientation matrix
  fRotMat (kNSD,kNSD),

  // crystal Cauchy stress
  fs_ij (kNSD),

  // crystal consistent tangent
  fc_ijkl (dSymMatrixT::NumValues(kNSD)),

  // anisotropic (cubic) part of elastic matrix
  fCanisoLat (dSymMatrixT::NumValues(kNSD)),   // lattice axes
  fCanisoBar (dSymMatrixT::NumValues(kNSD)),   // sample axes (Bbar)

  // Schmidt tensors in sample coords
  fZ (fNumSlip),
  fP (fNumSlip),

  // increm shearing rate
  fDGamma_n (fNumSlip),
  fDGamma   (fNumSlip),
  fdgam_save(fNumSlip),

  // resolve shear stress
  fTau (fNumSlip),

  // 2nd order identity tensor
  fISym (kNSD),

  // tensors in polar decomp
  fEigs (kNSD),
  fRe   (kNSD,kNSD),
  fUe   (kNSD),

  // work spaces
  fmatx1    (kNSD,kNSD),
  fmatx2    (kNSD,kNSD),
  fmatx3    (kNSD,kNSD),
  fmatx4    (kNSD,kNSD),
  fRank4    (dSymMatrixT::NumValues(kNSD)),
  fsymmatx1 (kNSD),
  fsymmatx2 (kNSD),

  // temp arrays
  fLHS   (fNumSlip),
  fRHS   (fNumSlip),
  fA     (fNumSlip),
  fB     (fNumSlip),
  farray (fNumSlip),

  // parameter in mid-point rule
  fTheta (1.0),

  // average stress
  fAvgStress (kNSD)
{
  // allocate additional space for Schmidt tensors
  for (int i = 0; i < fNumSlip; i++)
    {
      fZ[i].Allocate(3,3);
      fP[i].Allocate(3);
    }

  // allocate additional space for temp arrays
  for (int i = 0; i < fNumSlip; i++)
    {
      fA[i].Allocate(kNSD);
      fB[i].Allocate(kNSD);
      farray[i].Allocate(kNSD,kNSD);
    }

  // set anisotropic (cubic) part of elasticity matrix (lattice axes)
  fCanisoLat = 0.;
  fCanisoLat(0,0) = fCanisoLat(1,1) = fCanisoLat(2,2) = 1.;

  // set 2nd order unit tensor (sym matrix)
  fISym.Identity();
}

LocalCrystalPlast::~LocalCrystalPlast() {} 

int LocalCrystalPlast::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                       // fRotMat (const)
  d_size += kNSD*kNSD;                       // fFpi_n    (t_n)
  d_size += fNumSlip;                        // fDGamma_n (t_n)
  d_size += kNSD*kNSD;                       // fFpi      (t)
  d_size += fNumSlip;                        // fDGamma   (t)
  d_size += kNSD*kNSD;                       // fFe       (t)
  d_size += dim;                             // fs_ij     (t)
  d_size += fHardening->NumberOfVariables(); // Hard Vars at t_n and t

  // total # crystal variables per element (at all IP's)
  d_size *= NumIP() * fNumGrain;

  // averaged (aggregate) stress and moduli (at all IP's)
  d_size += NumIP() * dim;                   // fsavg_ij   (t)
  d_size += NumIP() * dim * dim;             // fcavg_ijkl (t)

  return d_size;
}

const dSymMatrixT& LocalCrystalPlast::s_ij()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // compute state, stress and moduli 
  if (fStatus == GlobalT::kFormRHS)
    {
      // reset iteration counter to check NLCSolver and state convergence
      if (CurrIP() == 0) 
	{
	  fIterCount = 0;
	  fIterState = 0;
	}

      // initialize average stress and moduli at IP
      fsavg_ij = 0.0;
      fcavg_ijkl = 0.0;

      // total deformation gradient
      // fFtot_n = fContinuumElement.FEManager().LastDeformationGradient();
      // fFtot = fContinuumElement.FEManager().DeformationGradient();
      fFtot_n = DeformationGradient(fLocLastDisp); 
      fFtot = DeformationGradient(fLocDisp);

      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // Schmidt tensors in sample coordinates (Bbar configuration)
	  for (int i = 0; i < fNumSlip; i++)
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	      fP[i].Symmetrize(fZ[i]);
	    }
	  
	  // fCanisoLat -> fCanisoBar (for anisotropic elasticity)
	  if (!fElasticity->IsIsotropic()) FFFFC_3D(fCanisoBar, fCanisoLat, fRotMat);

	  // compute crystal state
	  SolveCrystalState();
	  
	  // compute crystal Cauchy stress
	  CrystalS_ij();
	  
	  // compute plastic deformation gradient
	  FPInverse();

	  // compute crystal moduli
	  fc_ijkl = 0.;
          CrystalC_ijkl_Elastic();    // elastic part
	  CrystalC_ijkl_Plastic();    // plastic part

	  // add stress and moduli to corresponding averaged quantities
	  fsavg_ij.AddScaled(1./fNumGrain, fs_ij);
	  fcavg_ijkl.AddScaled(1./fNumGrain, fc_ijkl);
	}
    }

  // return averaged stress
  return fsavg_ij;
}

const dMatrixT& LocalCrystalPlast::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // return averaged moduli
  return fcavg_ijkl;
}

void LocalCrystalPlast::FormRHS(const dArrayT& dgamma, dArrayT& rhs)
{
  // incremental plastic deformation gradient ^ (-1)
  DeltaFPInverse(dgamma);

  // elastic right Cauchy-Green tensor 
  fCeBar.MultQTBQ(fDFpi, fCeBarTr);

  // resolved shear stress
  ResolveShearStress(fCeBar);

  // compute residual
  for (int i = 0; i < fNumSlip; i++)
    {
      // rhs[i] = fTau[i] - fKinetics->Psi(dgamma[i]/fdt, i); 
      rhs[i] = dgamma[i] - fdt * fKinetics->Phi(fTau[i], i);
    }
}

// form Jacobian for local Newton iteration
void LocalCrystalPlast::FormLHS(const dArrayT& dgamma, dMatrixT& lhs)
{
  // incremental plastic deformation gradient ^ (-1)
  //  DeltaFPInverse(dgamma);

  // elastic right Cauchy-Green tensor 
  //  fCeBar.MultQTBQ(fDFpi, fCeBarTr);

  // preliminary computations
  fdeltaI = 0.5 * fMatProp[1] * (fCeBar.Trace() - 3.) - fMatProp[0] + fMatProp[2];
  fCeBarTr.ToMatrix(fmatx1);
  fmatx4.MultATB(fDFpi, fmatx1); 
  dDFpidDGamma(dgamma, farray);

  // initialize jacobian
  lhs = 0.;

  // contribution of dTau/dDGamma to lhs(i,j)
  for (int i = 0; i < fNumSlip; i++)
    {
      // dTau/dCe
      dTaudCe(fZ[i], fP[i], fsymmatx1); 

      for (int j = 0; j < fNumSlip; j++)
	{
	  // dCe/dDGamma
	  dMatrixT& matx = farray[j];
	  fmatx3.MultAB(fmatx4, matx);
	  fsymmatx2.Symmetrize(fmatx3);
	  fsymmatx2 *= 2.;

	  // dTau/dCe : dCe/dDGamma
	  fsymmatx1.ToMatrix(fmatx1);
	  fsymmatx2.ToMatrix(fmatx2);
	  lhs(i,j) = dArrayT::Dot(fmatx1, fmatx2);
	}
    }

  // contribution of dPsi/dGamma to lhs(i,i) 
  /*for (int i = 0; i < fNumSlip; i++)
    {
      lhs(i,i) -= (fKinetics->DPsiDGamdot(dgamma[i]/fdt, i) / fdt);
      }*/

  for (int i = 0; i < fNumSlip; i++)
    {
      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      for (int j = 0; j < fNumSlip; j++) lhs(i,j) *= -tmp;
      lhs(i,i) += 1.;
    }

  // contribution of dRes/dHard*dHard/dDGamma to lhs(i,j)
  //  if (fAlgorCode == 3)
  //    {
  //      const dArrayT& dHdDGam = fHardening->ComputeHardQnts();
  //      for (int i = 0; i < fNumSlip; i++)
  //	{
  //	  double dResdHi = -fdt * fKinetics->DPhiDIso(fTau[i], i);
  //	  for (int j = 0; j < fNumSlip; j++) lhs(i,j) += dResdHi * dHdDGam[j];
  //	}
  //    }
}

void LocalCrystalPlast::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // update state at each integration point and ...
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // ... at each crystal
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // update state
	  fFpi_n      = fFpi;
	  fDGamma_n   = fDGamma;
          fHardening->UpdateHistory();
	}
    }
}

void LocalCrystalPlast::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // reset state at each integration point and ...
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // ... at each crystal
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // reset state
	  fFpi      = fFpi_n;
	  fDGamma   = fDGamma_n;
          fHardening->ResetHistory();
	}
    }
}

int LocalCrystalPlast::NumOutputVariables() const {return kNumOutput;}

void LocalCrystalPlast::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Allocate(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void LocalCrystalPlast::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int elem  = CurrElementNumber();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // Von Mises stress of the aggregate
  output[0] = sqrt(fsymmatx1.Deviatoric(fsavg_ij).ScalarProduct())/sqrt23;
  // cerr << " S_eq = " << output[0] << endl;
  //output[0] /= (fHardening->MaterialProperties())[1];

  // compute averaged equivalent stress
  if (elem == 0 && intpt == 0) fAvgStress = 0.0;
  fAvgStress.AddScaled(1./(NumIP()*NumElements()), fsavg_ij);
  // cout << " elem = " << elem << "   intpt = " << intpt << endl;
  // cout << "    fsavg_ij = " << endl << fsavg_ij << endl;
  // cout << "    fAvgStress = " << endl << fAvgStress << endl;
  if (elem == (NumElements()-1) && intpt == (NumIP()-1))
     cerr << " step # " << fContinuumElement.FEManager().StepNumber() 
          << "    S_eq_avg = " 
          << sqrt(fsymmatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23/4.0 << endl; 

  // iteration counter for nlcsolver and state
  output[1] = fIterCount;
  output[2] = fIterState;

  // compute texture of aggregate, if requested
  const int& step = fContinuumElement.FEManager().StepNumber();
  const int& nsteps = fContinuumElement.FEManager().NumberOfSteps();

  if (fmod(step, fODFOutInc) == 0 || step == nsteps)
    {
      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // fetch crystal data 
	  LoadCrystalData(element, intpt, igrn);
      
	  // texture: rotation tensor from fFe
	  PolarDecomp();

	  // texture: compute new crystal orientation matrix
	  fmatx1.MultAB(fRe, fRotMat);

	  // texture: compute Euler angles (in radians)
	  dArrayT& angles = fangles[igrn];
	  fLatticeOrient->RotMatrixToAngles(fmatx1, angles);
	}

      // write texture at IP/ELE
      fLatticeOrient->WriteTexture(elem, intpt, fNumGrain, step, fangles);
    }
}

void LocalCrystalPlast::Print(ostream& out) const
{
  // inherited
  PolyCrystalMatT::Print(out);

  // print slip kinetics data
  out << "    Kinetics of crystal slip (local crystal plast.)\n";
  out << "       Kinetics law. . . . . . . . . . . . . . . = " << fKinEqnCode << "\n";
  fKinetics->Print(out);

  // print slip hardening data
  out << "    Crystal slip hardening (local crystal plast.)\n";
  out << "       Hardening law . . . . . . . . . . . . . . = " << fHardCode << "\n";
  fHardening->Print(out);
}

void LocalCrystalPlast::PrintName(ostream& out) const
{
  // inherited
  PolyCrystalMatT::PrintName(out);

  // output model name
  out << "    Local crystal plasticity equations\n";

  // output detailed name of the model
  fSlipGeometry->PrintName(out);
  fElasticity->PrintName(out);
  fKinetics->PrintName(out);
  fHardening->PrintName(out);
}

GlobalT::SystemTypeT LocalCrystalPlast::TangentType() const
{
  return GlobalT::kNonSymmetric;
}

const dArrayT& LocalCrystalPlast::GetResolvedShearStress() const 
{ return fTau; }

const dArrayT& LocalCrystalPlast::GetIncrSlipShearStrain() const 
{ return fDGamma; }

/* PROTECTED MEMBER FUNCTIONS */

void LocalCrystalPlast::InitializeCrystalVariables()
{
  // initialize state at each element and ...
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);

      // ... at each integration point and ...
      for (int intpt = 0; intpt < NumIP(); intpt++)
	{
	  // load aggregate data at integration point
	  LoadAggregateData(element, intpt);

	  // initialilize average stress and moduli
	  fsavg_ij = 0.;
	  fcavg_ijkl = 0.;

	  // ... at each crystal
	  for (int igrn = 0; igrn < fNumGrain; igrn++)
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

	      // shear rates on slip systems
	      fDGamma_n = 0.;
	      fDGamma   = 0.;

	      // elastic deformation gradient
	      fFe.Identity();

	      // crystal Cauchy stress
	      fs_ij = 0.;

	      // hardening variables
	      fHardening->InitializeHardVariables(); 
	    }
	}
    }
}

void LocalCrystalPlast::LoadCrystalData(ElementCardT& element,
					int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int stressdim   = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 4*kNSD*kNSD + stressdim + 2*fNumSlip + fHardening->NumberOfVariables();
  int blockPerAgg = stressdim + stressdim*stressdim;
  int dex = intpt*(fNumGrain*blockPerGrn + blockPerAgg) + igrain*blockPerGrn;

  fRotMat.Set    (kNSD,kNSD, &d_array[dex             ]);   
  fFpi_n.Set     (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fDGamma_n.Set  (fNumSlip,  &d_array[dex += kNSD*kNSD]);
  fFpi.Set       (kNSD,kNSD, &d_array[dex += fNumSlip ]);
  fDGamma.Set    (fNumSlip,  &d_array[dex += kNSD*kNSD]);
  fFe.Set        (kNSD,kNSD, &d_array[dex += fNumSlip ]);     
  fs_ij.Set      (kNSD,      &d_array[dex += kNSD*kNSD]);
  fHardening->LoadHardData(stressdim, dex, d_array);
}

void LocalCrystalPlast::LoadAggregateData(ElementCardT& element, int intpt)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int stressdim = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 4*kNSD*kNSD + stressdim + 2*fNumSlip + fHardening->NumberOfVariables();
  int blockPerAgg = stressdim + stressdim*stressdim;
  int dex = intpt*(fNumGrain*blockPerGrn + blockPerAgg) + fNumGrain*blockPerGrn;

  fsavg_ij.Set   (kNSD,                &d_array[dex             ]);   
  fcavg_ijkl.Set (stressdim,stressdim, &d_array[dex += stressdim]);     
}

void LocalCrystalPlast::SetSlipKinetics()
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

    case SlipKinetics::kHaasen:          // Haasen power law (iso)
      fKinetics = new HaasenKinetics(*this);
      break;

    default:
      throwRunTimeError("LocalCrystalPlast::SetSlipKinetics: Bad fKinEqnCode");
    }
  if (!fKinetics) throwMemoryError("LocalCrystalPlast::SetSlipKinetics");
}

void LocalCrystalPlast::SetSlipHardening()
{
  // read hardening code model
  fInput >> fHardCode;

  // select slip hardening law
  switch(fHardCode)
    {
    case SlipHardening::kHard_L1:           // Voce's model (iso)
      fHardening = new VoceHardening(*this);
      break;
      
    case SlipHardening::kHard_L2:           // latent type hard law
      //fHardening = new LatentHardening(*this);
      throwRunTimeError("LocalCrystalPlast::SetSlipHardening: Not implemented");
      break;

    case SlipHardening::kHard_L3:           // Haasen's model (iso)
      fHardening = new HaasenHardening(*this);
      break;

    default:
      throwRunTimeError("LocalCrystalPlast::SetSlipHardening: Bad fHardCode");
    }
  if (!fHardening) throwMemoryError("LocalCrystalPlast::SetSlipHardening");
}

void LocalCrystalPlast::SolveCrystalState()
{
  // flag to track convergence of crystal state
  bool stateConverged = false;

  // counters for subincrementation (continuation method)
  int subIncr = 1;
  int totSubIncrs = 1;

  for(;;)
    {
      // time step
      fdt = fContinuumElement.FEManager().TimeStep() * (float)subIncr / (float)totSubIncrs;

      // deformation gradients
      fmatx1.SetToCombination(1., fFtot, -1., fFtot_n);
      //fFt_n.SetToCombination(1., fFtot_n, ((float)subIncr-1.0)/(float)totSubIncrs, fmatx1);
      fFt.SetToCombination(1., fFtot_n, (float)subIncr/(float)totSubIncrs, fmatx1);

      // save current solution  
      if (subIncr > 1 && stateConverged) {
         fdgam_save = fDGamma;
         fHardening->SaveCurrentSolution();
      }

      // iterate to compute crystal state
      IterateOnCrystalState(stateConverged, subIncr);

      // check convergence; if "stateConverged=false", use continuation method
      if (!stateConverged) 
	{
	  subIncr = 2 * subIncr - 1;
	  totSubIncrs = 2 * totSubIncrs;
          if (totSubIncrs > pow(2, 30)) 
              throwRunTimeError("LocalCrystalPlast::SolveCrystalState: totSubIncrs > 2^30");
          if (subIncr > 1) {
             fDGamma = fdgam_save;
             fHardening->RestoreSavedSolution();
          }
	}
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
      else
	break;
    }
}

void LocalCrystalPlast::IterateOnCrystalState(bool& stateConverged, int subIncr)
{
  // trial deformation quantities
  TrialDeformation();

  // initial guess for first sub-increment
  if (subIncr == 1)
    {
      // forward gradient estimate for DGamma
      ForwardGradientEstimate();

      // explicit estimate for hardness
      fHardening->ExplicitUpdateHard();
    }
  
  // iterate for crystal state
  int iterState = 0;
  int ierr = 0;

  int fAlgorCode = 2;
  switch(fAlgorCode)
    {
    case 1:
      while (++iterState <= fMaxIterState && !stateConverged)
	{
	  try
	    {
	      // iter level 1: solve for fDGamma; Hardness = constant
	      SolveForDGamma(ierr);
	      if (ierr != 0) {
		writeWarning("LocalCrystalPlast::SolveCrystalState:");
		writeWarning("   Will use continuation method");
		return;
	      }
	      
	      // iter level 2: solve for Hardness; fDGamma = constant
	      fHardening->ImplicitUpdateHard();
	      
	      // check convergence of state
	      stateConverged = (fHardening->Converged(fTolerState));
	    }

	  catch(int code)
	    {
	      break;
	    }
	}
      break;

    case 2:
      while (++iterState <= fMaxIterState && !stateConverged)
	{
          try
	    {
	      // iter level 1: solve for fDGamma; Hardness = constant
	      SolveForDGamma(ierr);
	      if (ierr != 0) {
		// writeWarning("LocalCrystalPlast::SolveCrystalState:
		//   ierr!=0 in SolveForDGamma: Will try continuation method");
		return;
	      }
	      
	      // iter level 2: solve for Hardness; fDGamma = constant
	      fHardening->ImplicitSolveHard();
	      
	      // check convergence of state
	      stateConverged = (Converged(fTolerState) && fHardening->Converged(fTolerState));
	    }
	  
          catch(int code)
	    {
              // writeWarning("LocalCrystalPlast::SolveCrystalState: 
              //    exception thrown & caugth at SolveForDGamma, Will try continuation method");
	      break;
	    }
	}
      break;

    case 3: // NEED TO BE IMPROVED !!!
      while (++iterState <= fMaxIterState && !stateConverged)
	{
	  // inner loop: Newton solve for Hardness; fDGamma = constant
	  fHardening->ImplicitSolveHard();

	  // outer loop: compute new Newton iterate for fDGamma
	  fSolver->SolveOneIteration(fSolverPtr, fDGamma, ierr);
	  if (ierr == 0) stateConverged = true;
	  if (ierr == 1) continue;
	  if (ierr == 2) break;
	}
      break;

    default:
      throwRunTimeError("LocalCrystalPlast::SolveCrytalState: Bad fAlgorCode");
    }

  // check if did not converge in max iterations
  if (!stateConverged && iterState > fMaxIterState) {
    // writeWarning("LocalCrystalPlast::SolveCrystalState: 
    //   didn't converge in maxIters, Will try continuation method");
    return;
  }
  
  // update iteration counter for state
  fIterState = max(fIterState, --iterState);
}

void LocalCrystalPlast::TrialDeformation()
{
  // trial elastic deformation gradient
  fFeTr.MultAB(fFt, fFpi_n);

  // trial elastic right Cauchy-Green tensor
  fCeBarTr.MultATA(fFeTr);
}

/*void LocalCrystalPlast::ForwardGradientEstimate()
{
  // fDGamma_ini(t) computed for incr = 1
  // fDGamma_ini(t) = FDGamma(t_n) for incr > 1
  if (ftime <= fdt) 
    {
      // shear stress on slip systems
      ResolveShearStress(fCeBarTr);
      
      // initial shear deformation
      for (int i = 0; i < fNumSlip; i++)
	{
	  //fDGamma[i] = 0.0;
	  fDGamma[i] = fdt * fKinetics->Phi(fTau[i], i);
	}
    }
}*/

void LocalCrystalPlast::ForwardGradientEstimate()
{
  // reference to hardening/kinetics material properties
  const dArrayT& propH = fHardening->MaterialProperties();
  const dArrayT& propKE = fKinetics->MaterialProperties();
  double m = propKE[0];

  // some local tensors
  dMatrixT fFe_n (kNSD);
  dSymMatrixT fbe_n (kNSD);
  dSymMatrixT fsigma_n (kNSD);

  // elasticity tensors at t_n
  // fFtot_n = DeformationGradient(fLocLastDisp); 
  // fFtot_n = fContinuumElement.FEManager().LastDeformationGradient();
  fFe_n.MultAB(fFtot_n, fFpi_n);
  fbe_n.MultAAT(fFe_n);
  
  // rate of deformation tensor at B_n
  fmatx1.Inverse(fFtot_n);
  fmatx2.MultAB(fFt, fmatx1);                // F_rel
  fsymmatx1.MultATA(fmatx2);                 // D
  fsymmatx1.PlusIdentity(-1.);
  fsymmatx1 *= 0.5 / fdt;

  // Schmidt tensors at t_n
  fmatx1.Inverse(fFe_n);
  for (int i = 0; i < fNumSlip; i++)
    {
      fmatx2.MultAB(fFe_n, fZ[i]);
      farray[i].MultAB(fmatx2, fmatx1);      // Z_n
      fA[i].Symmetrize(farray[i]);           // P_n
    }

  // resolved shear stress and Kirchhoff stress at t_n
  fsymmatx2.MultATA(fFe_n);
  ResolveShearStress(fsymmatx2);
  fsigma_n.MultQBQT(fFe_n, fSBar);

  // sym matrices (D, be_n) -> full matrix representation
  fsymmatx1.ToMatrix(fmatx1);
  fbe_n.ToMatrix(fmatx2);

  for (int i = 0; i < fNumSlip; i++) 
    {
      fsymmatx1.MultQBQT(fmatx2, fA[i]);
      fsymmatx1 *= 2.*fMatProp[0];

      fsigma_n.ToMatrix(fmatx4);
      fmatx3.MultAB(farray[i], fmatx4);
      fsymmatx2.Symmetrize(fmatx3);

      fA[i].ToMatrix(fmatx3);
      fsymmatx1.AddCombination(fMatProp[1]*dArrayT::Dot(fmatx2, fmatx3), fbe_n,
			       2., fsymmatx2);
      fsymmatx1.ToMatrix(fmatx3);

      //if (ftime <= fdt)
      if (ftime >= 0.0)
	fDGamma[i] = dArrayT::Dot(fmatx3, fmatx1) * fdt;
      else
	fDGamma[i] = fDGamma_n[i] * (1. + fdt / (m*fTau[i]) * dArrayT::Dot(fmatx3, fmatx1));

      for (int j = 0; j < fNumSlip; j++)
	{
	  fA[j].ToMatrix(fmatx4);
	  // if (ftime <= fdt) {
	  if (ftime >= 0.0) {
	    fLHS(i, j) = dArrayT::Dot(fmatx3, fmatx4);
	    // if (i == j) fLHS(i, j) += propH[0];
	    if (i == j) fLHS(i, j) += fHardening->HardeningModulus();
	  }
	  else {
	    fLHS(i, j) = fDGamma_n[i] / (m*fTau[i]) * dArrayT::Dot(fmatx3, fmatx4);
	    if (i == j) fLHS(i, j) += 1.;
	  }
	}

    }

  // solve for initial estimate of dgamma
  fLHS.LinearSolve(fDGamma);

  // norm of initial estimate of dgamma
  fMagDGam0 = fDGamma.Magnitude();
} 

bool LocalCrystalPlast::Converged(double toler)
{
  // check convergence on dgamma
  bool test = ( fabs(fMagDGam-fMagDGam0) < toler*fMagDGam0 );

  // if did not converge, reset norm
  if (!test) fMagDGam0 = fMagDGam;

  return test;
}

void LocalCrystalPlast::CrystalS_ij()
{
  // incremental plastic deformation gradient ^ (-1)
  DeltaFPInverse(fDGamma);

  // Elastic Deformation Gradient
  fFe.MultAB(fFeTr, fDFpi);

  // elastic rigth Cauchy-Green tensor
  fCeBar.MultATA(fFe);

  // 2nd Piola Kirchhoff Stress
  CrystalPKIIStress(fCeBar);

  // Cauchy Stress
  fs_ij.MultQBQT(fFe, fSBar);
  fs_ij /= fFe.Det();
}

void LocalCrystalPlast::FPInverse()
{
  // Incremental Plastic Deformation Gradient 
  DeltaFPInverse(fDGamma);

  // Current Plastic Deformation Gradient
  fFpi.MultAB(fFpi_n, fDFpi);

  // Normalize Plastic Deformation Gradient
  if (fFpi.Det() <= 0.0) 
    throwRunTimeError("LocalCrystalPlast::FPInverse: det(Fpi) < 0");
  fFpi /= pow(fFpi.Det(), 1./3.);
}

void LocalCrystalPlast::CrystalC_ijkl_Elastic()
{
  // elastic left Cauchy-Green tensor (b_e)
  fsymmatx1.MultAAT(fFe);

  // I_b tensor
  Set_I_b_Tensor(fsymmatx1, fc_ijkl);
  
  // b_e (x) b_e
  fRank4.Outer(fsymmatx1, fsymmatx1);

  // isotropic elastic contribution to spatial moduli
  fc_ijkl *= 2.*fMatProp[0];
  fc_ijkl.AddScaled(fMatProp[1], fRank4);

  // anisotropic (cubic) elastic contribution to spatial moduli
  if (!fElasticity->IsIsotropic()) 
    {
      FFFFC_3D(fRank4, fCanisoBar, fFe);
      fc_ijkl.AddScaled(-2.*fMatProp[2], fRank4);
    }
}

/* implements a rough approximation to consistent tangent */
/*void LocalCrystalPlast::CrystalC_ijkl_Plastic()
{
  // temp matrices
  dMatrixT fcp (dSymMatrixT::NumValues(kNSD));
  dMatrixT fcei (dSymMatrixT::NumValues(kNSD));

  // inverse of elasti part
  fcei = MatrixInversion(fc_ijkl);

  // compute plastic part
  fcp = 0.;
  fmatx1.Inverse(fFe);
  for (int i = 0; i < fNumSlip; i++)
    {
      fmatx2.MultAB(fFe, fZ[i]);
      fmatx3.MultAB(fmatx2, fmatx1);     // Z = Fe Z0 Fe^(-1)
      fsymmatx1.Symmetrize(fmatx3);      // P = (Z)_sym
      fRank4.Outer(fsymmatx1, fsymmatx1);
      double tmp = fdt * fKinetics->dPhidTau(fTau[i], i);
      fcp.AddScaled(tmp, fRank4);         // dDp/dSig
    }

  // tangent stiffness
  fRank4.SetToCombination(1., fcei, fdt, fcp);
  fc_ijkl = MatrixInversion(fRank4);
  }*/

/*void LocalCrystalPlast::CrystalC_ijkl_Plastic()
{
  // inverse of plastic right C-G tensor
  dSymMatrixT fCpi(kNSD);
  fCpi.MultAAT(fDFpi);

  // shear stress on slip systems
  ResolveShearStress(fCeBar);
  
  // SECOND TERM : (Fe_tr) (dDGamma/dCe_tr) (Fe_tr)^T
  FormLHS(fDGamma, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  for (int i = 0; i < fNumSlip; i++)
    {
      dTaudCe(fZ[i], fP[i], fsymmatx1);
      fA[i].MultQBQT(fDFpi, fsymmatx1);

      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      fA[i] *= -tmp;
    }

  for (int i = 0; i < fNumSlip; i++) {
    fB[i] = 0.;
    for (int j = 0; j < fNumSlip; j++)
      fB[i].AddScaled(-Jaci(i,j), fA[j]);
  }

  for (int i = 0; i < fNumSlip; i++) 
    {
      dSymMatrixT& matx = fB[i]; 
      fsymmatx2.MultQBQT(fFeTr, matx);
      fB[i] = fsymmatx2;
    }

  // FIRST TERM : (Fe_tr) (dS/dDGamma) (Fe_tr)^T
  fCeBarTr.ToMatrix(fmatx2);
  fCpi.ToMatrix(fmatx3);
  for (int i = 0; i < fNumSlip; i++)
    {
      dMatrixT& matx = farray[i];      // dDFpi/dDGamma from FormLHS()
      fmatx1.MultABT(matx, fDFpi);
      fsymmatx1.Symmetrize(fmatx1);
      fsymmatx1 *= 2.;

      fsymmatx1.ToMatrix(fmatx1);
      fmatx4.MultAB(fmatx1, fmatx2);
      fmatx1.MultAB(fmatx4, fmatx3);
      fsymmatx2.Symmetrize(fmatx1);
      fsymmatx2 *= 2.*fMatProp[0];
      
      fsymmatx1.ToMatrix(fmatx1);
      fsymmatx2.AddCombination(fdeltaI, fsymmatx1, 
                               0.5*fMatProp[1]*dArrayT::Dot(fmatx2,fmatx1), fCpi);
      fA[i].MultQBQT(fFeTr, fsymmatx2);
    }
	
  // PLASTIC CONTRIBUTION TO MODULI : 2 * Sum(fA[i] (x) fB[i])
  for (int i = 0; i < fNumSlip; i++)
    {
      fRank4.Outer(fA[i], fB[i]);
      fc_ijkl.AddScaled(2., fRank4);
    }
}*/

void LocalCrystalPlast::CrystalC_ijkl_Plastic()
{
  // compute inverse of fDFpi
  fDFp.Inverse(fDFpi);

  // shear stress on slip systems
  ResolveShearStress(fCeBar);
  
  // SECOND TERM : (Fe_tr) (dDGamma/dCe_tr) (Fe_tr)^T
  FormLHS(fDGamma, fLHS);
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

void LocalCrystalPlast::SolveForDGamma(int& ierr)
{
  // solve for incremental shear strain
  fSolver->Solve(fSolverPtr, fDGamma, ierr);

  if (ierr != 0) {
    writeWarning("LocalCrystalPlast::SolveForDGamma:");
    writeWarning("   Convergence problems in NLCSolver");
    return;
  }

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());

  // norm of DGamma
  fMagDGam = fDGamma.Magnitude();
}

void LocalCrystalPlast::ResolveShearStress(const dSymMatrixT& Ce)
{
  // 2nd Piola Kirchhoff Stress
  CrystalPKIIStress(Ce);

  // Resolve Shear Stress on Slip Systems
  Ce.ToMatrix(fmatx2);
  fSBar.ToMatrix(fmatx3);
  fmatx1.MultAB(fmatx2, fmatx3); 
  for(int i = 0; i < fNumSlip; i++)
    {
      fTau[i] = dArrayT::Dot(fmatx1, fZ[i]);
    }
}

/* mid-point rule to compute fDFpi */
void LocalCrystalPlast::DeltaFPInverse(const dArrayT& dgamma)
{
  // Lambda: Sum_(i=1)^(i=fNumslip) dgamma[i]*Z[i]
  fmatx1 = 0.;
  for (int i = 0; i < fNumSlip; i++)
    fmatx1.AddScaled(dgamma[i], fZ[i]);

  // I-theta*Lambda
  fmatx2.Identity();
  fmatx2.AddScaled(-fTheta, fmatx1);

  // (I+(1-theta)*Lambda)^(-1)
  fmatx3.Identity();
  fmatx3.AddScaled((1.-fTheta), fmatx1);
  fmatx3.Inverse();

  // incremental plastic deformation gradient ^ (-1)
  fDFpi.MultAB(fmatx3, fmatx2);
}

void LocalCrystalPlast::CrystalPKIIStress(const dSymMatrixT& Ce)
{
  // elastic Green strain 
  fEeBar.SetToCombination(0.5, Ce, -0.5, fISym);

  // isotropic contribution to 2nd P-K stress: 2*mu*Ee + lambda*tr(Ee)*I
  fSBar.SetToCombination(2.*fMatProp[0], fEeBar, fMatProp[1]*fEeBar.Trace(), fISym);

  // anisotropic contribution to 2nd P-K stress: -2*beta*Caniso[Ee]
  if (!fElasticity->IsIsotropic()) 
    {
      dSymMatrixT symmatx3(kNSD), symmatx4(kNSD);
      //symmatx3.DoubleOffDiags(fEeBar);
      //fCanisoBar.Multx(symmatx3, symmatx4);
      symmatx4.A_ijkl_B_kl(fCanisoBar, fEeBar);
      fSBar.AddScaled(-2.*fMatProp[2], symmatx4 );
    }
}

/* dDFpi/dDGamma for mid-point rule */
void LocalCrystalPlast::dDFpidDGamma(const dArrayT& dgamma, ArrayT<dMatrixT>& array)
{
  // Lambda: Sum_(i=1)^(i=fNumslip) dgamma[i]*Z[i]
  fmatx1 = 0.;
  for (int i = 0; i < fNumSlip; i++)
    fmatx1.AddScaled(dgamma[i], fZ[i]);
  
  // [I+(1-theta)*Lambda]^(-1)
  fmatx3.Identity();
  fmatx3.AddScaled((1.-fTheta), fmatx1);
  fmatx3.Inverse();

  // [theta*I+(1-theta)*fDFpi]
  fmatx2.Identity(fTheta);
  fmatx2.AddScaled((1.-fTheta), fDFpi);

  // dDFpi/dDGamma
  for (int i = 0; i < fNumSlip; i++)
    {
      fmatx1.MultAB(fmatx3, fZ[i]);
      array[i].MultAB(fmatx1, fmatx2);
      array[i] *= -1.;
    }
}

void LocalCrystalPlast::dTaudCe(const dMatrixT& Z, const dSymMatrixT& P, 
				dSymMatrixT& symmatx)
{
  // isotropic contribution to dTau/dCe
  P.ToMatrix(fmatx2);
  fCeBar.ToMatrix(fmatx3);
  fmatx1.MultAB(fmatx2, fmatx3);
  symmatx.Symmetrize(fmatx1);
  symmatx *= 2. * fMatProp[0];
  symmatx.AddCombination(fdeltaI, P, 0.5*fMatProp[1]*dArrayT::Dot(fmatx3,fmatx2), fISym);

  // anisotropic contribution to dTau/dCe
  if (!fElasticity->IsIsotropic())
    {
      dSymMatrixT symmatx3(kNSD), symmatx4(kNSD), symmatx5(kNSD);
      //symmatx3.DoubleOffDiags(fCeBar);
      //fCanisoBar.Multx(symmatx3, symmatx4);
      symmatx4.A_ijkl_B_kl(fCanisoBar, fCeBar);
      symmatx4.ToMatrix(fmatx1);
      fmatx2.MultAB(Z, fmatx1);
      symmatx3.Symmetrize(fmatx2);             // (X*C[CeBar])_sym
    
      fmatx1.MultAB(fmatx3, Z);
      //symmatx4.Symmetrize(fmatx1);
      //symmatx5.DoubleOffDiags(symmatx4);
      //fCanisoBar.Multx(symmatx5, symmatx4);   // C[(CeBar*X)_sym]
      symmatx5.Symmetrize(fmatx1);
      symmatx4.A_ijkl_B_kl(fCanisoBar, symmatx5);
      
      symmatx.AddCombination(-fMatProp[2], symmatx3, -fMatProp[2], symmatx4);
    }
}

void LocalCrystalPlast::PolarDecomp()
{
  // elastic right Cauchy-Green tensor
  fCeBar.MultATA(fFe);

  // principal values of right stretch tensor
  fCeBar.PrincipalValues(fEigs);
  for (int i = 0; i < kNSD; i++)
    {
      if (fEigs[i] < 0.) throwRunTimeError("LocalCrystalPlast::PolarDecomp: eigs < 0");
      fEigs[i] = sqrt(fEigs[i]);
    }

  // temporaries
  double i1 = fEigs[0] + fEigs[1] + fEigs[2];
  double i2 = fEigs[0]*fEigs[1] + fEigs[0]*fEigs[2] + fEigs[1]*fEigs[2];
  double i3 = fEigs[0]*fEigs[1]*fEigs[2];
  double D = i1*i2 - i3;
  if (D < 0.0) throwRunTimeError("LocalCrystalPlast::PolarDecomp: D < 0");
  
  // elastic right stretch tensor
  fsymmatx1.MultAB(fCeBar, fCeBar);
  fUe.SetToCombination(-1., fsymmatx1, (i1*i1-i2), fCeBar, i1*i3, fISym);
  fUe /= D;

  // inverse of fUe
  fsymmatx1.SetToCombination(1., fCeBar, -i1, fUe, i2, fISym);
  fsymmatx1 /= i3;
  fsymmatx1.ToMatrix(fmatx1);

  // rotation tensor
  fRe.MultAB(fFe, fmatx1); 
}

const dMatrixT& LocalCrystalPlast::DeformationGradient(const LocalArrayT& disp)
{ 
  return F(disp); 
}

void LocalCrystalPlast::Set_I_b_Tensor(const dSymMatrixT& b, dMatrixT& c)
{
  double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
  double z13, z14, z15, z16, z17, z18, z19, z20, z21;
  
  z1 = b[0];
  z2 = b[1];
  z3 = b[2];
  z4 = b[3];
  z5 = b[4];
  z6 = b[5];
  z7 = z1*z1;
  z8 = z1*z2;
  z9 = z2*z2;
  z10 = z1*z3;
  z11 = z2*z3;
  z12 = z3*z3;
  z13 = z1*z4;
  z14 = z2*z4;
  z15 = z3*z4;
  z16 = z4*z4;
  z17 = z1*z5;
  z18 = z2*z5;
  z19 = z3*z5;
  z20 = z4*z5;
  z21 = z5*z5;
  z1 = z1*z6;
  z2 = z2*z6;
  z3 = z3*z6;
  z4 = z4*z6;
  z5 = z5*z6;
  z6 = z6*z6;
  z11 = z11 + z16;
  z10 = z10 + z21;
  z3 = z20 + z3;
  z18 = z18 + z4;
  z13 = z13 + z5;
  z8 = z6 + z8;
  z11 = 0.5*z11;
  z10 = 0.5*z10;
  z3 = 0.5*z3;
  z18 = 0.5*z18;
  z13 = 0.5*z13;
  z8 = 0.5*z8;
  
  //{{z7, z6, z21, z5, z17, z1}, 
  // {z6, z9, z16, z14, z4, z2}, 
  // {z21, z16, z12, z15, z19, z20}, 
  // {z5, z14, z15, z11, z3, z18}, 
  // {z17, z4, z19, z3, z10, z13}, 
  // {z1, z2, z20, z18, z13, z8}}
  
  c(0,0) = z7;
  c(0,1) = c(1,0) = z6;
  c(0,2) = c(2,0) = z21;
  c(0,3) = c(3,0) = z5;
  c(0,4) = c(4,0) = z17;
  c(0,5) = c(5,0) = z1;
  c(1,1) = z9;
  c(1,2) = c(2,1) = z16;
  c(1,3) = c(3,1) = z14;
  c(1,4) = c(4,1) = z4;
  c(1,5) = c(5,1) = z2;
  c(2,2) = z12;
  c(2,3) = c(3,2) = z15;
  c(2,4) = c(4,2) = z19;
  c(2,5) = c(5,2) = z20;
  c(3,3) = z11;
  c(3,4) = c(4,3) = z3;
  c(3,5) = c(5,3) = z18;
  c(4,4) = z10;
  c(4,5) = c(5,4) = z13;
  c(5,5) = z8;
}

void LocalCrystalPlast::FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F)
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
