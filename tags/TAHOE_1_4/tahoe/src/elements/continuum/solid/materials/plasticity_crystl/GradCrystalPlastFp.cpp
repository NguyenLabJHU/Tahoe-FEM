/* $Id: GradCrystalPlastFp.cpp,v 1.16 2004-04-14 19:52:20 ebmarin Exp $ */
#include "GradCrystalPlastFp.h"
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

/* useful constant */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

/* to debug */
const bool XTAL_MESSAGES = false;
const int ELprnt = 0;

GradCrystalPlastFp::GradCrystalPlastFp(ifstreamT& in, const FSMatSupportT& support) :
  LocalCrystalPlastFp(in, support),  
  fLocInitX   (ContinuumElement().InitialCoordinates()),
  fLocCurrX   (LocalArrayT::kCurrCoords),
//  fLocInitXIP (NumIP(), NumSD()),
//  fLocCurrXIP (NumIP(), NumSD()),
  fGradTool   (new GradientTools_C(NumIP(), NumSD())),
  fFpIP       (NumIP()),    
  fFpC        (kNSD),
  fGradFp     (kNSD),
  fCurlFpT    (kNSD),
  fKe_n       (kNSD,kNSD),  
  fKe         (kNSD,kNSD),  
  fXe         (kNSD,kNSD),  
  fnormFp0    (NumIP()),
  fnormHard0  (NumIP()),
  fnormFp     (NumIP()),
  fnormHard   (NumIP()),
  fMatx4      (kNSD,kNSD),
  fX_IP       (NumSD())
{
  // check number of grains
  if (fNumGrain != 1) 
    throwRunTimeError("GradCrystalPlastFp::GradCrystalPlastFp: NumGrain != 1");

  // may need to check element node number(?)
  int nodes = fLocInitX.NumberOfNodes();

  // check number of IPs used (if #sd=2 -> #ip=4; if #sd=3 -> #ip=8)
  if (NumSD() == 2)
    if (NumIP() != 4) 
      throwRunTimeError("GradCrystalPlastFp::GradCrystalPlastFp: NumSD=2 && NumIP!=4");
  else if (NumSD() == 3)
    if (NumIP() != 8) 
      throwRunTimeError("GradCrystalPlastFp::GradCrystalPlastFp: NumSD=3 && NumIP!=8");
  else
      throwRunTimeError("GradCrystalPlastFp::GradCrystalPlastFp: NumSD!=2 or 3");

  // allocate space for ...
  // ... Fp values at integration points
  for (int i = 0; i < NumIP(); i++)
    fFpIP[i].Dimension(kNSD,kNSD);

  // ... spatial gradients of Fp (note: kNSD instead of NumSD())
  for (int i = 0; i < kNSD; i++)
    fGradFp[i].Dimension(kNSD,kNSD);

  // ... current coordinates
  fLocCurrX.Dimension(fLocInitX.NumberOfNodes(), NumSD());

  // ... ip coordinates
  fLocInitXIP.Dimension(NumIP(), NumSD());
  fLocCurrXIP.Dimension(NumIP(), NumSD());
}

GradCrystalPlastFp::~GradCrystalPlastFp() {} 

int GradCrystalPlastFp::NumVariablesPerElement()
{
  int d_size = 0;
  int strdim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                        // fRotMat   (const)
  d_size += kNSD*kNSD;                        // fFp_n     (t_n)
  d_size += kNSD*kNSD;                        // fFe_n     (t_n)
  d_size += kNSD*kNSD;                        // fKe_n     (t_n)
  d_size += kNSD*kNSD;                        // fFp       (t)
  d_size += kNSD*kNSD;                        // fFe       (t)
  d_size += kNSD*kNSD;                        // fKe       (t)
  d_size += kNSD*kNSD;                        // fXe       (t)
  d_size += fNumSlip;                         // fDGamma   (t)
  d_size += strdim;                           // fs_ij     (t)
  d_size += strdim*strdim;                    // fc_ijkl   (t)
  d_size += strdim;                           // fC        (t) (bookeeping)
  d_size += strdim*strdim;                    // fcBar_ijkl(t) (bookeeping)
  d_size += fHardening->NumberOfVariables();  // hard vars at t and t_n

  // total # variables per element
  d_size *= fNumGrain * NumIP();

  return d_size;
}

const dSymMatrixT& GradCrystalPlastFp::s_ij()
{
  // fetch current element
  ElementCardT & element = CurrentElement();

  // only one grain per integration point
  int igrn = 0;

  // time step
  fdt = fFSMatSupport.TimeStep();

  // compute crystal stresses (all IPs at once - elastic predictor at first iter)
  if (fFSMatSupport.RunState() == GlobalT::kFormRHS && CurrIP() == 0)
    {
       if (fFSMatSupport.IterationNumber() <= -1)
         {
           for (int intpt = 0; intpt < NumIP(); intpt++)
             {
               // deformation gradient
               Compute_Ftot_3D(fFtot, intpt);

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

               // elastic tensors
               fFpi.Inverse(fFp_n);
               fFe.MultAB(fFtot, fFpi);
               fCeBar.MultATA(fFe);

               // 2nd Piola Kirchhoff Stress
               fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
               fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

               // Cauchy Stress
               fs_ij.MultQBQT(fFe, fSBar);
               fs_ij /= fFe.Det();
             }
         }
       else
         {
           // state at all integration points
           SolveCrystalState();

           // crystal stresses 
           for (int intpt = 0; intpt < NumIP(); intpt++)
	     {
               // deformation gradient at intpt
               Compute_Ftot_3D(fFtot, intpt);

	       // fetch crystal data
	       LoadCrystalData(element, intpt, igrn);
	  
	       // compute crystal Cauchy stress
	       CrystalS_ij();
 	     }
         }
    }

  // fetch crystal data
  LoadCrystalData(element, CurrIP(), igrn);

  // return crystal stress
  return fs_ij;
}

const dMatrixT& GradCrystalPlastFp::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // only one grain per integration point
  int igrn = 0;

  // recover local data
  LoadCrystalData(element, intpt, igrn);

  // elasticity matrix in Bbar configuration
  if (!fElasticity->IsIsotropic())
     {
        fElasticity->ComputeModuli(fRank4);
        FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
     }
  else
        fElasticity->ComputeModuli(fcBar_ijkl);

  if (fFSMatSupport.IterationNumber() <= 0)
    {
      // elastic crystal stiffness
      FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
    }
  else
    {
      // Schmidt tensor in Bbar configuration (sample axes)
      for (int i = 0; i < fNumSlip; i++)
         fZ[i].MultQBQT(fRotMat, fZc[i]);

      // compute crystal moduli (consistent tangent)
      CrystalC_ijkl();
    }

  // return crystal moduli
  return fc_ijkl;
}

void GradCrystalPlastFp::UpdateHistory()
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
      fFp_n     = fFp;
      fFe_n     = fFe;
      fKe_n     = fKe;
      //fDGamma_n = fDGamma;
      fHardening->UpdateHistory();
    }
}

void GradCrystalPlastFp::ResetHistory()
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
      fFp     = fFp_n;
      fFe     = fFe_n;
      fKe     = fKe_n;
      //fDGamma = fDGamma_n;
      fHardening->ResetHistory();
    }
}

int GradCrystalPlastFp::NumOutputVariables() const {return kNumOutput;}

void GradCrystalPlastFp::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void GradCrystalPlastFp::ComputeOutput(dArrayT& output)
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
  output[0] = sqrt(fSymMatx1.Deviatoric(fs_ij).ScalarProduct())/sqrt23;
  // cerr << " S_eq = " << output[0] << endl;

  // compute averaged equivalent stress
  if (elem == 0 && intpt == 0) fAvgStress = 0.0;
  fAvgStress.AddScaled(1./(NumIP()*NumElements()), fs_ij);
  if (elem == (NumElements()-1) && intpt == (NumIP()-1))
     cerr << " step # " << fFSMatSupport.StepNumber()
          << "    S_eq_avg = " 
          << sqrt(fSymMatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23
          << "    Savg_12 = " << fAvgStress(0,1) << endl; 

  // iteration counters
  output[1] = fIterCount;
  output[2] = fIterState;

  // compute euler angles
  int step = fFSMatSupport.StepNumber();
  int nsteps = fFSMatSupport.NumberOfSteps();

  if (fmod(double(step), fODFOutInc) == 0 || step == nsteps)
  {
    // texture: rotation tensor from fFe
    fSpecD.PolarDecomp(fFe, fRe, fUe, false);
    //PolarDecomp();
    
    // texture: compute new crystal orientation matrix
    fMatx1.MultAB(fRe, fRotMat);

    // texture: compute Euler angles (in radians)
    dArrayT& angles = fangles[igrn];
    fLatticeOrient->RotMatrixToAngles(fMatx1, angles);

    // write euler angles at IP/ELE
    fLatticeOrient->WriteTexture(group, elem, intpt, fNumGrain, step, fangles);
  }
}

void GradCrystalPlastFp::Print(ostream& out) const
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

void GradCrystalPlastFp::PrintName(ostream& out) const
{
  // inherited
  PolyCrystalMatT::PrintName(out);

  // output model name
  out << "    gradient crystal plasticity equations - Fp\n";
  fSlipGeometry->PrintName(out);
  fKinetics->PrintName(out);
  fHardening->PrintName(out);
}

/* PROTECTED MEMBER FUNCTIONS */

void GradCrystalPlastFp::SetSlipKinetics()
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
      throwRunTimeError("GradCrystalPlastFp::SetSlipKinetics: Bad fKinEqnCode");
    }
  if (!fKinetics) throwMemoryError("GradCrystalPlastFp::SetSlipKinetics");
}

void GradCrystalPlastFp::SetSlipHardening()
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
      throwRunTimeError("GradCrystalPlastFp::SetSlipHardening: Bad fHardCode");
    }
  if (!fHardening) throwMemoryError("GradCrystalPlastFp::SetSlipHardening");
}

void GradCrystalPlastFp::InitializeCrystalVariables(ElementCardT& element)
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
	  
	  // plastic/elastic deformation gradients
	  fFp_n.Identity();
	  fFe_n.Identity();
	  fFp.Identity();
	  fFe.Identity();

	  // elastic curvatures
	  fKe_n = 0.;
	  fKe   = 0.;
	  
	  // stress associated with lattice curvature
	  fXe = 0.;
	      
          // shear rates on slip systems
	  fDGamma = 0.;
	      
	  // crystal Cauchy stress and moduli
	  fs_ij = 0.;
	  fc_ijkl = 0.;

	  // hardening variables
	  fHardening->InitializeHardVariables();
	}
}

void GradCrystalPlastFp::LoadCrystalData(ElementCardT& element,
				       int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int sdim        = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 8*kNSD*kNSD + fNumSlip + 2*sdim + 2*sdim*sdim
                    + fHardening->NumberOfVariables();
  int dex         = intpt*fNumGrain*blockPerGrn + igrain*blockPerGrn;

  fRotMat.Set    (kNSD,kNSD, &d_array[dex             ]);   
  fFp_n.Set      (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fFe_n.Set      (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);   
  fKe_n.Set      (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);   
  fFp.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fFe.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fKe.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);   
  fXe.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);   
  fDGamma.Set    (fNumSlip,  &d_array[dex += kNSD*kNSD]);
  fs_ij.Set      (kNSD,      &d_array[dex += fNumSlip ]);
  fc_ijkl.Set    (sdim,sdim, &d_array[dex += sdim     ]);
  fC.Set         (kNSD,      &d_array[dex += sdim*sdim]);   
  fcBar_ijkl.Set (sdim,sdim, &d_array[dex += sdim     ]);
  fHardening->LoadHardData(sdim*sdim, dex, d_array);
}

void GradCrystalPlastFp::SolveCrystalState()
{
  // fetch current element
  ElementCardT& element = CurrentElement();

  // reset iteration counter to check NLCSolver
  fIterState = 0;
  fIterCount = 0;

  // only one grain per integration point
  int igrn = 0;

  // get coordinates of IPs at initial configuration
  for (int intpt = 0; intpt < NumIP(); intpt++)
    { 
//      dummy = fFSMatSupport.Interpolate(fLocInitX, fX_IP, intpt);
      ContinuumElement().IP_Interpolate(fLocInitX, fX_IP, intpt);
      for (int j = 0; j < NumSD(); j++) fLocInitXIP(intpt,j) = fX_IP[j];
    }

  // shape function derivarives dNa/dXinit at center (to compute fGradFp)
  fGradTool->ComputeGDNa(fLocInitXIP);

  // elastic curvature
  for (int intpt = 0; intpt < NumIP(); intpt++) 
    {
      LoadCrystalData(element, intpt, igrn);
      fFpIP[intpt] = fFp_n;
    }
  LatticeCurvature(element, igrn);

  // get reference to hardening material properties
  const dArrayT& prop = fHardening->MaterialProperties();

  dArrayT array(3);

  for (int intpt = 0; intpt < NumIP(); intpt++) 
    {
      // deformation gradient at integration point IP
      Compute_Ftot_3D(fFt, intpt);
      Compute_Ftot_last_3D(fFtot_n, intpt);

      // fetch crystal data
      LoadCrystalData(element, intpt, igrn);
      
      // Schmidt tensor in sample coordinates
      for (int i = 0; i < fNumSlip; i++) 
         {
            fZ[i].MultQBQT(fRotMat, fZc[i]);
            fRotMat.Multx(fSlipMc[i], fSlipM[i]);
         }

      // elasticity matrix in Bbar configuration
      if (!fElasticity->IsIsotropic())
         {
            fElasticity->ComputeModuli(fRank4);
            FFFFC_3D(fcBar_ijkl, fRank4, fRotMat);
         }
      else
            fElasticity->ComputeModuli(fcBar_ijkl);

      // right Cauchy-Green tensor
      fC.MultATA(fFt);

      // conjugate stress to elastic curvature & incompatibility stress-like qnt
      //fXe.SetToScaled(prop[4]*prop[5]*fMatProp[0], fKe_n); 
      fXe.SetToScaled(prop[4]*prop[5]*fMatProp[0], fKe); 
      for (int i = 0; i < fNumSlip; i++)
         {
            // Achyara's suggestion
            fKe.MultTx(fSlipM[i], array);
            fHardening->fTauInc[i] = prop[6]*fMatProp[0]*array.Magnitude();

            // Gurtin's suggestion
            //fKe_n.Multx(fSlipM[i], array);
            //fHardening->fTauInc[i] = prop[6]*fMatProp[0]*dArrayT::Dot(fSlipM[i], array);
         }
      
      // initial guess for Fp
      InitialEstimateForFp();

      // initial guess for hardness
      InitialEstimateForHardening();

      // initial norms for Fp and hardness
      fnormFp0[intpt]   = fFpNorm0;
      fnormHard0[intpt] = fHardening->Magnitude();

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << " IP # (Estimates) " << intpt << "\n"
              << "    fFtot: \n" << fFt << endl;
         cout << "    fFp  : \n" << fFp << endl;
         cout << "    fKe_n: \n" << fKe_n << endl;
         cout << "    fXe  : \n" << fXe << endl;
         cout << "    fTauKin   : \n" << fHardening->fTauKin << endl;
         cout << "    fDGamma   : \n" << fDGamma << endl;
         cout << "    normFp0   : " << fnormFp0[intpt] << endl;
         cout << "    normHard0 : " << fnormHard0[intpt] << endl;
      }
    }

  // iterate for the state in the element (all IP's at once)
  int ierr;
  int iterState;
  bool stateConverged = false;
  for (iterState = 0; iterState < fMaxIterState; iterState++)
    {
      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << "\n\n***Iterate for state" << endl;
         cout << "   iterState # " << iterState << "\n\n";
         cout << "   === Compute Fp ===" << endl;
      }

      // compute Fp at each integration point
      for (int intpt = 0; intpt < NumIP(); intpt++) 
	{
	  // fetch crystal data
	  LoadCrystalData(element, intpt, igrn);

	  // Schmidt tensors and shear back stress
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
	    }
	  
	  // solve for Fp
	  SolveForPlasticDefGradient(ierr);
	  
	  // bookeeping to compute lattice curvature
	  //fFpIP[intpt] = fFp;

	  // norm for Fp
          fnormFp[intpt] = fFpNorm;

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             cout << " IP # " << intpt << "\n";
             cout << "    No iters : " << fIterCount << "\n";
             cout << "    fFp      : \n" << fFp << endl;
             cout << "    fTauKin  : \n" << fHardening->fTauKin << endl;
             cout << "    fDGamma  : \n" << fDGamma << endl;
             cout << "    normFp   : " << fnormFp[intpt] << endl;
          }
	}
      
      // elastic curvature
      //LatticeCurvature(element, igrn);

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt)
         cout << "\n\n   === Compute Hardness ===" << endl;

      // compute hardness
      for (int intpt = 0; intpt < NumIP(); intpt++) 
	{
	  // fetch crystal data
	  LoadCrystalData(element, intpt, igrn);
	  
	  // Schmidt tensor in sample coordinates
	  for (int i = 0; i < fNumSlip; i++) 
	    {
	      fZ[i].MultQBQT(fRotMat, fZc[i]);
              fRotMat.Multx(fSlipMc[i], fSlipM[i]);
	    }

          // conjugate stress to elastic curvature & incompatibility stress-like qnt
          fXe.SetToScaled(prop[4]*prop[5]*fMatProp[0], fKe); 
          for (int i = 0; i < fNumSlip; i++)
             {
                // Achyara's suggestion
                fKe.MultTx(fSlipM[i], array);
                fHardening->fTauInc[i] = prop[6]*fMatProp[0]*array.Magnitude();

                // Gurtin's suggestion
                //fKe.Multx(fSlipM[i], array);
                //fHardening->fTauInc[i] = prop[6]*fMatProp[0]*dArrayT::Dot(fSlipM[i], array);
             }
      
          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             cout << " IP # " << intpt << "\n";
             cout << "    fXe      : \n" << fXe << endl;
          }

	  // compute hardness
	  SolveForHardening();

	  // norm for hardness
	  fnormHard[intpt] = fHardening->Magnitude();

          if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
             //cout << " IP # " << intpt << "\n";
             //cout << "    fXe      : \n" << fXe << endl;
             cout << "    fTauKin  : \n" << fHardening->fTauKin << endl;
             cout << "    normHard : " << fnormHard[intpt] << endl;
          }
	}

      // check state convergence
      stateConverged = CheckConvergenceOfState();

      if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
         cout << "\n\n   === Convergence Check ===\n";
         cout << " stateConverged: " << stateConverged << endl;
         for (int i=0; i<NumIP(); i++)
           cout << " IP # " << i  << endl
                << "     normHard0 : " << fnormHard0[i] << endl
                << "     normHard  : " << fnormHard[i] << endl
                << "     normFp0   : " << fnormFp0[i] << endl
                << "     normFp    : " << fnormFp[i] << endl
                << "     | diff normFp   | = " << fabs(fnormFp[i]-fnormFp0[i]) << endl
                << "     | diff normHard | = " << fabs(fnormHard[i]-fnormHard0[i]) << endl;
      }

      if (stateConverged) break;

      // reset norms
      fnormFp0   = fnormFp;
      fnormHard0 = fnormHard;
    }


  // check if did not converge in max iterations
  if (!stateConverged) {
    writeWarning("... in GradCrystalPlastFp::SolveCrystalState: iters > MaxIters\n ...... will throw 'eBadJacobianDet' to cut dtime");
    throw ExceptionT::kBadJacobianDet;
  }

  // update iteration count for state
  if (stateConverged) fIterState = max(fIterState, iterState);
}

void GradCrystalPlastFp::SolveForPlasticDefGradient(int& ierr)
{
  // 9x1 array form of Fp
  Rank2ToArray9x1(fFp, fFpArray);

  // uses "kind of" continuation method based on rate sensitivity exponent
  fKinetics->SetUpRateSensitivity();

  ierr = 0;
  do {
       // current value for rate sensitivity exponent
       fKinetics->ComputeRateSensitivity();
 
       // solve for incremental shear strain
       try { fSolver->Solve(fSolverPtr, fFpArray, ierr); }
       catch(ExceptionT::CodeT code)
           {
             fKinetics->RestoreRateSensitivity();
             writeWarning("... in GradCrystalPlastFp::SolveForPlasticDefGradient: caugth exception\n ...... will throw 'eBadJacobianDet' to cut dtime");
             throw ExceptionT::kBadJacobianDet;
           }

       // throw ExceptionT::xception if problems in NCLSolver
       if (ierr != 0) {
          fKinetics->RestoreRateSensitivity();
          writeWarning("... in GradCrystalPlastFp::SolveForPlasticDefGradient: ierr!=0\n ...... will throw 'eBadJacobianDet' to cut dtime");
          throw ExceptionT::kBadJacobianDet;
       }

     } while (!fKinetics->IsMaxRateSensitivity());

  // update iteration count from NLCSolver
  fIterCount = max(fIterCount, fSolver->GetIterationCount());

  // recover 3x3 matrix form of Fp
  Rank2FromArray9x1(fFp, fFpArray);
  fFp /= pow(fFp.Det(), 1./3.);

  // norm of Fp
  fFpNorm = sqrt(fFp.ScalarProduct());
}

void GradCrystalPlastFp::ResolveShearStress()
{
  // 2nd Piola Kirchhoff stress
  fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
  fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

  // resolve shear stress on slip systems
  fCeBar.ToMatrix(fMatxCe);
  fSBar.ToMatrix(fMatxSb);
  fMatx1.MultAB(fMatxCe, fMatxSb); 
  for (int i = 0; i < fNumSlip; i++)
     fTau[i] = dArrayT::Dot(fMatx1, fZ[i]);

  // resolve back stress on slip systems
  fMatx1.MultAB(fMatxCe, fXe);
  //fMatx1.SetToScaled(1., fXe);
  for (int i = 0; i < fNumSlip; i++) 
     fHardening->fTauKin[i] = dArrayT::Dot(fMatx1, fZ[i]);
}

/* PRIVATE MEMBER FUNCTIONS */

bool GradCrystalPlastFp::CheckConvergenceOfState() const
{
  int i = 0;
  bool test = true;
  while (i < NumIP() && test)
    {
      test = ( fabs(fnormFp[i]-fnormFp0[i]) <= fTolerState*fnormFp0[i] &&
	       fabs(fnormHard[i]-fnormHard0[i]) <= fTolerState*fnormHard0[i] );
      i++;
    }
  return test;
}

void GradCrystalPlastFp::LoadCrystalCurvature(ElementCardT& element, 
					    int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int sdim        = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 8*kNSD*kNSD + fNumSlip + 2*sdim + 2*sdim*sdim
                    + fHardening->NumberOfVariables();
  int dex         = intpt*fNumGrain*blockPerGrn + igrain*blockPerGrn;
  int blockToKe   = 6*kNSD*kNSD;

  // current elastic curvature
  fKe.Set (kNSD,kNSD, &d_array[dex += blockToKe]);
}

// dislocation tensor aBar^p
void GradCrystalPlastFp::LatticeCurvature(ElementCardT& element, int igrn)
{
  // interpolate Fp from IPs to center of IP-elemen
  fGradTool->InterpolateTensorAtCenter(fFpIP, fFpC); 

  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
     cout << "\n\n   === Lattice Curvature ===" << endl;
        cout << " FpC " << endl
             << fFpC << "    Det(FpC) = " << fFpC.Det() << endl;
  }

  // normalize Fp at center such that Det(Fp)=1
  fFpC /= pow(fFpC.Det(), 1./3.);

  // Curl of Fp^T at center of IP-element
  fGradTool->CurlOfTensorAtCenter(fFpIP, fCurlFpT);

  // curvature tensor Ke at center of IP-element (dislocation tensor aBar^p)
  fMatx1.MultAB(fFpC, fCurlFpT);
  fMatx2.SetToScaled(1./fFpC.Det(), fMatx1);

  if (XTAL_MESSAGES && CurrElementNumber() == ELprnt) {
     cout << " Center " << "\n";
     cout << "    fCurlFpT: \n" << fCurlFpT << endl;
     cout << "    Curvature tensor fKe: \n" << fMatx2 << endl;
   }

  // assigned Ke_center to all IPs
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // fetch crystal curvature
      LoadCrystalCurvature(element, intpt, igrn);

      // curvature tensor fKe at IP (dislocation tensor aBar^p)
      fKe = fMatx2;
    }
}

void GradCrystalPlastFp::AddGradTermToLHS(dMatrixT& lhs, const dMatrixT& matx)
{
  // term:
  //  -SUM <dGamDot/dKin*Z (x) {[Ce*Xe*Z^T + trace(Ce*Xe*Z^T)]*matx^T}>
  //  where: matx = fFpi (LocalJac);  matx = fIMatx (moduli)
  fMatx1.MultAB(fMatxCe, fXe);
  for (int i = 0; i < fNumSlip; i++)
    {
      fMatx2.MultABT(fMatx1, fZ[i]);
      fMatx3.SetToCombination(1., fMatx2, fMatx2.Trace(), fIMatx);
      //fMatx4.MultABT(fMatx3, fFpi);
      fMatx4.MultABT(fMatx3, matx);
      OuterProduct9x9(fZ[i], fMatx4, fRankIV_1);
      lhs.AddScaled(-fKinetics->DPhiDKin(fTau[i],i), fRankIV_1);  
    }
}

void GradCrystalPlastFp::AddGradTermToC_ijkl()
{
  // term : 4 * SUM [dGamDot/dKin * A (x) (Z*Xe^T)_s]
  // where: A = (Y*SBar)_s + 0.5*(Cijkl*(CeBar*Y)_s)
  for (int i = 0; i < fNumSlip; i++)
    {
      fMatx1.MultABT(fZ[i], fXe);
      fSymMatx1.Symmetrize(fMatx1);
      fRank4.Outer(fA[i], fSymMatx1);
      fcBar_ijkl.AddScaled(4.*fKinetics->DPhiDKin(fTau[i],i), fRank4);
    }
}
