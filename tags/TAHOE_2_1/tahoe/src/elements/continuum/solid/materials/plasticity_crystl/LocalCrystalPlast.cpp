/* $Id: LocalCrystalPlast.cpp,v 1.21 2004-07-15 08:29:07 paklein Exp $ */
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

#include "Utils.h"
#include "SpectralDecompT.h"
#include "ContinuumElementT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 3;
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

/* to debug */
const bool XTAL_MESSAGES = false;
const int IPprnt = 1;

LocalCrystalPlast::LocalCrystalPlast(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("local_crystal_plasticity"),
//	PolyCrystalMatT(in, support),  

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

  // crystal consistent tangent (Bbar config)
  fcBar_ijkl (dSymMatrixT::NumValues(kNSD)),

  // Schmidt tensors in sample coords
  fP (fNumSlip),

  // increm shearing rate
  fDGamma_n (fNumSlip),
  fdgam_save(fNumSlip),

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
  fsymmatx3 (kNSD),

  // temp arrays
  fLHS   (fNumSlip),
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
      fP[i].Dimension(3);
    }

  // allocate additional space for temp arrays
  for (int i = 0; i < fNumSlip; i++)
    {
      fA[i].Dimension(kNSD);
      fB[i].Dimension(kNSD);
      farray[i].Dimension(kNSD,kNSD);
    }

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

int LocalCrystalPlast::NumberOfUnknowns() const { return fNumSlip; }

const dSymMatrixT& LocalCrystalPlast::s_ij()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // compute state, stress and moduli 
  if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
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
      Compute_Ftot_last_3D(fFtot_n);
      Compute_Ftot_3D(fFtot);

      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
	  // recover local data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // Schmidt tensor in sample coordinates (Bbar configuration)
	  for (int i = 0; i < fNumSlip; i++)
            {
              fZ[i].MultQBQT(fRotMat, fZc[i]);
              fP[i].Symmetrize(fZ[i]);
            }

          // elasticity matrix in Bbar configuration
          if (!fElasticity->IsIsotropic())
             {
                fElasticity->ComputeModuli(fRank4);
                FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
             }
          else
                fElasticity->ComputeModuli(fcBar_ijkl);
          //cout << " fcBar_ijkl " << endl << fcBar_ijkl << endl;

          if (XTAL_MESSAGES && CurrIP() == IPprnt) {
             cout << " Step # " << 
                      fFSMatSupport->StepNumber() 
                  << "    Iter # " << 
                      fFSMatSupport->IterationNumber() 
                  << endl;
          }

          // compute crystal Cauchy stress and consistent tangent
          // global iterations start at iter = -1
          if (fFSMatSupport->StepNumber() >= 0 &&
              fFSMatSupport->IterationNumber() <= -1)
             {
               // deformation gradient
               fmatx1.SetToCombination(1., fFtot, -1., fFtot_n);
               fFt.SetToCombination(1., fFtot_n, 100.0 / 100.0, fmatx1);

               // elastic tensors
               fFe.MultAB(fFt, fFpi_n);
               fCeBar.MultATA(fFe);
 
               // 2nd Piola Kirchhoff Stress at t_n
               fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
               fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

               // Cauchy Stress
               fs_ij.MultQBQT(fFe, fSBar);
               fs_ij /= fFe.Det();

               // elastic crystal stiffness
               FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
             }
          else 
             {
	       // compute crystal state
	       SolveCrystalState();
	  
       	       // compute plastic deformation gradient
	       FPInverse();

	       // compute crystal Cauchy stress
	       CrystalS_ij();
	  
	       // compute crystal moduli (consistent tangent)
	       CrystalC_ijkl();  
             }

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
  ResolveShearStress();

  // compute residual
  for (int i = 0; i < fNumSlip; i++)
    {
      //rhs[i] = fTau[i] - fKinetics->Psi(dgamma[i]/fdt, i); 
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

  // term: X_b = 2*DFp*dDFpi/dDgamma_b = 2*DFp*(-Z_b)
  // (assumes backward Euler integration for FpDot)
  fDFp.Inverse(fDFpi);     
  //for (int i = 0; i < fNumSlip; i++)
  //  {
  //    farray[i].MultAB(fDFp, fZ[i]);
  //    farray[i] *= -2.0;       
  //  }
  ArrayT<dMatrixT> array(fNumSlip);
  for (int i = 0; i < fNumSlip; i++)
    {
      array[i].Dimension(kNSD,kNSD);
    }
  dDFpidDGamma(dgamma, array);
  for (int i = 0; i < fNumSlip; i++)
    {
      farray[i].MultAB(fDFp, array[i]);
      farray[i] *= 2.0;       
    }

  // initialize local Jacobian
  lhs = 0.;

  // term: dTau/dDGamma = dTau/dCe : dCe/dDGamma
  fCeBar.ToMatrix(fmatx3);
  fSBar.ToMatrix(fmatx4);
  for (int i = 0; i < fNumSlip; i++)
    {
      // dTau/dCe
      fmatx1.MultAB(fmatx3, fZ[i]);   // (CeBar*Z)
      fmatx2.MultAB(fZ[i], fmatx4);   // (Z*SBar)
      fsymmatx1.Symmetrize(fmatx1);   // (CeBar*Z)_s
      fsymmatx2.Symmetrize(fmatx2);   // (Z*SBar)_s

      fsymmatx3.A_ijkl_B_kl(fcBar_ijkl, fsymmatx1); // (Cijkl*(CeBar*Z)_s)
      fsymmatx2.AddScaled(0.5, fsymmatx3);   // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
      fsymmatx2.ToMatrix(fmatx2);

      for (int j = 0; j < fNumSlip; j++)
	{
	  // dCe/dDGamma
	  dMatrixT& Xbeta = farray[j];   // Xb
	  fmatx1.MultAB(fmatx3, Xbeta);  // (CeBar*Xb)
	  fsymmatx1.Symmetrize(fmatx1);  // (CeBar*Xb)_s
	  fsymmatx1.ToMatrix(fmatx1);

	  // A_ab = dTau/dDGamma = (Z*SBar)_s+0.5*(Cijkl*(CeBar*Z)_s):(CeBar*Xb)_s
	  lhs(i,j) = dArrayT::Dot(fmatx2, fmatx1);
	}
    }

  // local Jacobian: SUM_b ( delta_ab - dt * dPhi_a/dTau_a * A_ab )
  for (int i = 0; i < fNumSlip; i++)
    {
      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      for (int j = 0; j < fNumSlip; j++) lhs(i,j) *= -tmp;
      lhs(i,i) += 1.;
    }

  /*for (int i = 0; i < fNumSlip; i++)
    {
      lhs(i,i) -= (fKinetics->DPsiDGamdot(dgamma[i]/fdt, i) / fdt);
    }*/

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
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void LocalCrystalPlast::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int group = ContinuumElement().ElementGroupNumber();
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
     cerr << " step # " << fFSMatSupport->StepNumber() 
          << "    group # " << group
          << "    S_eq_avg = " 
          << sqrt(fsymmatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23 << endl; 

  // iteration counter for nlcsolver and state
  output[1] = fIterCount;
  output[2] = fIterState;

  // compute texture of aggregate, if requested
  int step = fFSMatSupport->StepNumber();
  int nsteps = fFSMatSupport->NumberOfSteps();

  if (fmod(double(step), fODFOutInc) == 0 || step == nsteps)
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
      fLatticeOrient->WriteTexture(group, elem, intpt, fNumGrain, step, fangles);
    }
}

GlobalT::SystemTypeT LocalCrystalPlast::TangentType() const
{
  return GlobalT::kNonSymmetric;
}

/* PROTECTED MEMBER FUNCTIONS */

void LocalCrystalPlast::InitializeCrystalVariables(ElementCardT& element)
{
	// current element number
	int elem = CurrElementNumber();

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

void LocalCrystalPlast::IterateOnCrystalState(bool& stateConverged, int subIncr)
{
  // initialze flag to track convergence
  stateConverged = false;

  if (XTAL_MESSAGES && CurrIP() == IPprnt) 
     cout << " Starting solution: fDGamma = \n" << fDGamma << endl;

  // trial deformation quantities
  TrialDeformation();

  // initial guess for first sub-increment
  if (subIncr == 1) {
    
      // forward gradient estimate for DGamma
      ForwardGradientEstimate();

      if (XTAL_MESSAGES && CurrIP() == IPprnt)
         cout << " Fwd Grad Guess : fDGamma = \n" << fDGamma << endl;

      // explicit estimate for hardness
      fHardening->ExplicitUpdateHard();
    }
  
  // iterate for crystal state
  int iterState = 0;
  int ierr = 0;

  int fAlgorCode = 1;
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
		writeWarning("LocalCrystalPlast::IterateOnCrystalState: ierr != 0 in SolveForDGamma -> subincrementation");
		cout << " elem # " << CurrElementNumber()
	             << ";  IP # " << CurrIP() << endl;
		cout << " sunIncr # " << subIncr << endl;
		return;
	      }
	      
	      // iter level 2: solve for Hardness; fDGamma = constant
	      fHardening->ImplicitUpdateHard();
	      
	      // check convergence of state
	      stateConverged = (Converged(fTolerState) && fHardening->Converged(fTolerState));
	      //stateConverged = (fHardening->Converged(fTolerState));
	    }

	  catch(ExceptionT::CodeT code)
	    {
               if (XTAL_MESSAGES) {
                  writeWarning("LocalCrystalPlast::IterateOnCrystalState: exception caugth at SolveForDGamma -> subincrementation"); 
		  cout << " elem # " << CurrElementNumber()
    	               << ";  IP # " << CurrIP() << endl;
		  cout << " subIncr # " << subIncr << endl;
               }
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
	  
          catch(ExceptionT::CodeT code)
	    {
               if (XTAL_MESSAGES) {
                  writeWarning("LocalCrystalPlast::IterateOnCrystalState: exception caugth at SolveForDGamma -> subincrementation"); 
                  cout << " IP # " << CurrIP() << endl;
               }
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
      throwRunTimeError("LocalCrystalPlast::IterateOnCrystalState: Bad fAlgorCode");
    }

  // check if did not converge in max iterations
  if (!stateConverged && iterState > fMaxIterState) {
     writeWarning("LocalCrystalPlast::IterateOnCrystalState: didn't converge in maxIters -> subincrementation");
     cout << " elem # " << CurrElementNumber()
          << ";  IP # " << CurrIP() << endl;
     cout << " sunIncr # " << subIncr << endl;
     return;
  }
  
  // update iteration counter for state
  if (stateConverged) fIterState = max(fIterState, --iterState);
}

void LocalCrystalPlast::TrialDeformation()
{
  // trial elastic deformation gradient
  fFeTr.MultAB(fFt, fFpi_n);

  // trial elastic right Cauchy-Green tensor
  fCeBarTr.MultATA(fFeTr);
}

void LocalCrystalPlast::ForwardGradientEstimate()
{
  // reference to hardening/kinetics material properties
  const dArrayT& propH = fHardening->MaterialProperties();
  const dArrayT& propKE = fKinetics->MaterialProperties();
  double m = propKE[0];
  double time = fFSMatSupport->Time();

  // some local tensors
  dMatrixT fFe_n (kNSD);
  dSymMatrixT fCeDot (kNSD);

  // elasticity tensors at t_n
  fFe_n.MultAB(fFtot_n, fFpi_n);
  fCeBar.MultATA(fFe_n);
 
  // rate : 0.5*L_v^p(CeBar)
  fmatx1.Inverse(fFtot_n);
  fmatx2.MultAB(fFt, fmatx1);                // frel
  fmatx3.Inverse(fmatx2);                    // (frel)^(-1)
  fsymmatx1.MultATA(fmatx2);                 
  fsymmatx1.PlusIdentity(-1.);
  fsymmatx1 *= 0.5 / fdt;                    // Edot=0.5*(frel^T frel - 1)/dt
  //fsymmatx2.MultQTBQ(fmatx3, fsymmatx1);     // Edot=frel^(-T)*Edot*frel^(-1)
  //fCeDot.MultQTBQ(fFe_n, fsymmatx2);         // Cedot = 0.5*L_v^p(CeBar) = Fe^T*Edot*Fe
  fCeDot.MultQTBQ(fFe_n, fsymmatx1);         // Cedot = 0.5*L_v^p(CeBar) = Fe^T*Edot*Fe

  // resolved shear stress
  ResolveShearStress();

  for (int i = 0; i < fNumSlip; i++)
    {
      fmatx1.MultAB(fmatx3, fZ[i]);   // (CeBar*Z)
      fmatx2.MultAB(fZ[i], fmatx4);   // (Z*SBar)
      fB[i].Symmetrize(fmatx1);       // (CeBar*Z)_s
      fA[i].Symmetrize(fmatx2);       // (Z*SBar)_s
    }

  // Form LHS
  fLHS = 0.;

  fCeDot.ToMatrix(fmatx3);
  double theta = 0.5;
  for (int i = 0; i < fNumSlip; i++)
    {
      fsymmatx3.A_ijkl_B_kl(fcBar_ijkl, fB[i]); // (Cijkl*(CeBar*Z)_s)
      fA[i].AddScaled(0.5, fsymmatx3);   // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
      fA[i].ToMatrix(fmatx2);

      double tmp =  2. * fdt * dArrayT::Dot(fmatx2, fmatx3);
      if (time >= 0.0) 
      //if (ftime <= fdt)
        fDGamma[i] = tmp;
      else
	fDGamma[i] = fDGamma_n[i] + theta*fabs(fDGamma_n[i])/(m*fabs(fTau[i]))*tmp;

      for (int j = 0; j < fNumSlip; j++)
	{
           fsymmatx1.SetToScaled(2.0, fB[j]);    // 2*(CeBar*Z)_s
	   fsymmatx1.ToMatrix(fmatx1);

	   // (Z*SBar)_s+0.5*(Cijkl*(CeBar*Z)_s):(2*CeBar*Z)_s
	   fLHS(i,j) = dArrayT::Dot(fmatx2, fmatx1);

	   if (time >= 0.0) {
	   //if (ftime <= fdt) {
	      //if (i == j) fLHS(i, j) += fHardening->HardeningModulus();
	      if (i == j) fLHS(i, j) += 1.;
	   }
	   else {
	      fLHS(i, j) *= theta*fabs(fDGamma_n[i])/(m*fabs(fTau[i]));
	      if (i == j) fLHS(i, j) += 1.;
	   }
	}
    }

  // solve for initial estimate of dgamma
  fLHS.LinearSolve(fDGamma);

  // norm of initial estimate of dgamma
  fMagDGam0 = fDGamma.Magnitude();
} 

void LocalCrystalPlast::RestoreSavedSolution()
{
  fDGamma = fdgam_save;
  fHardening->RestoreSavedSolution();
}

void LocalCrystalPlast::SaveCurrentSolution()
{
  fdgam_save = fDGamma;
  fHardening->SaveCurrentSolution();
}

bool LocalCrystalPlast::Converged(double toler)
{
  // check convergence on dgamma
  bool test = ( fabs(fMagDGam-fMagDGam0) < toler*fMagDGam0 );

  // if did not converge, reset norm
  if (!test) fMagDGam0 = fMagDGam;

  return test;
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

void LocalCrystalPlast::CrystalS_ij()
{
  // incremental plastic deformation gradient ^ (-1)
  DeltaFPInverse(fDGamma);

  // Elastic Deformation Gradient
  fFe.MultAB(fFeTr, fDFpi);

  // elastic rigth Cauchy-Green tensor
  fCeBar.MultATA(fFe);

  // 2nd Piola Kirchhoff Stress
  fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
  fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

  // Cauchy Stress
  fs_ij.MultQBQT(fFe, fSBar);
  fs_ij /= fFe.Det();
}

/*void LocalCrystalPlast::CrystalC_ijkl()
{
  // elastic contribution to spatial moduli
  FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);

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
      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      fcp.AddScaled(tmp, fRank4);         // dDp/dSig
    }

  // tangent stiffness
  fRank4.SetToCombination(1., fcei, fdt, fcp);
  fc_ijkl = MatrixInversion(fRank4);
}*/

void LocalCrystalPlast::CrystalC_ijkl()
{
  // incremental plastic deformation gradient ^ (-1)
  DeltaFPInverse(fDGamma);

  // elastic right Cauchy-Green tensor 
  fCeBar.MultQTBQ(fDFpi, fCeBarTr);

  // resolved shear stress
  ResolveShearStress();

  // SECOND TERM
  FormLHS(fDGamma, fLHS);
  dMatrixT JacInv = MatrixInversion(fLHS);

  fCeBar.ToMatrix(fmatx3);
  fSBar.ToMatrix(fmatx4);
  for (int i = 0; i < fNumSlip; i++)
    {
      // dTau/dCe
      fmatx1.MultAB(fmatx3, fZ[i]);   // (CeBar*Z)
      fmatx2.MultAB(fZ[i], fmatx4);   // (Z*SBar)
      fsymmatx1.Symmetrize(fmatx1);   // (CeBar*Z)_s
      fsymmatx2.Symmetrize(fmatx2);   // (Z*SBar)_s

      fsymmatx3.A_ijkl_B_kl(fcBar_ijkl, fsymmatx1); // (Cijkl*(CeBar*Z)_s)
      fsymmatx2.AddScaled(0.5, fsymmatx3);   // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)

      double tmp = fdt * fKinetics->DPhiDTau(fTau[i], i);
      fA[i].SetToScaled(tmp, fsymmatx2);
    }

  for (int i = 0; i < fNumSlip; i++) 
    {
      fB[i] = 0.;
      for (int j = 0; j < fNumSlip; j++)
         fB[i].AddScaled(JacInv(i,j), fA[j]);
    }

  // FIRST TERM
  for (int i = 0; i < fNumSlip; i++)
    {
      // dTau/dCe
      dMatrixT& Xbeta = farray[i];    // farray from formLHS
      fmatx1.MultAB(fmatx3, Xbeta);   // (CeBar*X)
      fmatx2.MultAB(Xbeta, fmatx4);   // (X*SBar)
      fsymmatx1.Symmetrize(fmatx1);   // (CeBar*X)_s
      fsymmatx2.Symmetrize(fmatx2);   // (X*SBar)_s

      fsymmatx3.A_ijkl_B_kl(fcBar_ijkl, fsymmatx1); // (Cijkl*(CeBar*X)_s)
      fsymmatx2.AddScaled(0.5, fsymmatx3);   // (X*SBar)_s + 0.5*(Cijkl*(CeBar*X)_s)
      fA[i].SetToScaled(1.0, fsymmatx2);
    }

  // consistent tangent in Bbar:  CepBar = CeBar + 2 * SUM(fA[i] (x) fB[i])
  for (int i = 0; i < fNumSlip; i++)
    {
      fRank4.Outer(fA[i], fB[i]);
      fcBar_ijkl.AddScaled(2.0, fRank4);
    }

  // consistent tangent in current config: c_ijkl = Fe_iI*Fe_jJ*Fe_kK*Fe_lL*CepBar_IJKL
  FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
  //cout << " fc_ijkl E-P (New) " << endl << fc_ijkl << endl;
}

void LocalCrystalPlast::SolveForDGamma(int& ierr)
{
  // uses "kind of" continuation method based on rate sensitivity exponent
  fKinetics->SetUpRateSensitivity();

  do {
       // current value for rate sensitivity exponent
       fKinetics->ComputeRateSensitivity();
 
       // solve for incremental shear strain
       try { fSolver->Solve(fSolverPtr, fDGamma, ierr); }
       catch(ExceptionT::CodeT code)
           {
             fKinetics->RestoreRateSensitivity();
             throw;
           }

       // return if problems in NLCSolver
       if (ierr != 0) {
          fKinetics->RestoreRateSensitivity();
          return;
       }
     } while (!fKinetics->IsMaxRateSensitivity());

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());

  // norm of DGamma
  fMagDGam = fDGamma.Magnitude();
}

void LocalCrystalPlast::ResolveShearStress()
{
  // 2nd Piola Kirchhoff Stress
  fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
  fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

  // Resolve Shear Stress on Slip Systems
  fCeBar.ToMatrix(fmatx3);
  fSBar.ToMatrix(fmatx4);
  fmatx1.MultAB(fmatx3, fmatx4); 
  for (int i = 0; i < fNumSlip; i++)
     fTau[i] = dArrayT::Dot(fmatx1, fZ[i]);
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

void LocalCrystalPlast::PolarDecomp()
{
  // polar decomposition
  SpectralDecompT fSpecD(kNSD);
  fSpecD.PolarDecomp(fFe, fRe, fUe, false);

/*
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
*/  
}

/* to be deleted */
void LocalCrystalPlast::dTaudCe(const dMatrixT& Z, const dSymMatrixT& P,
                                dSymMatrixT& symmatx)
{ 
#pragma unused(Z)
#pragma unused(P)
#pragma unused(symmatx)
}

void LocalCrystalPlast::CrystalC_ijkl_Elastic() { }

void LocalCrystalPlast::CrystalC_ijkl_Plastic() { }
