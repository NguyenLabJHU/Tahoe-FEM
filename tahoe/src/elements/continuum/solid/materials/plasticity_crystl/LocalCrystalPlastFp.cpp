/* $Id: LocalCrystalPlastFp.cpp,v 1.19 2004-07-15 08:29:07 paklein Exp $ */
#include "LocalCrystalPlastFp.h"
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
#include "CommunicatorT.h"

#include "ContinuumElementT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 3;

/* useful constant */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {"VM_stress", "IterNewton", "IterState"};

/* to output some messages */
const bool XTAL_MESSAGES = true;
const int IPprnt = 1;

LocalCrystalPlastFp::LocalCrystalPlastFp(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("local_crystal_plasticity_Fp"),
//  PolyCrystalMatT(in, support),  

  // penalty parameter for detFp
  fPenalty (1.0e+0),

  // elastic deformation gradient
  fFe_n (kNSD,kNSD),
  fFe   (kNSD,kNSD),

  // plastic deformation gradients 
  fFp_n    (kNSD,kNSD),
  fFp      (kNSD,kNSD),
  fFpi     (kNSD,kNSD),
  fFp_save (kNSD,kNSD),

  // right Cauchy-Green tensor
  fC       (kNSD), 

  // symmetric tensors in interm config	
  fCeBar   (kNSD),
  fEeBar   (kNSD),
  fSBar    (kNSD),

  // crystal consistent tangent in Bbar configuration
  fcBar_ijkl (dSymMatrixT::NumValues(kNSD)),

  // tensors/objects in polar decomp
  fEigs  (kNSD),
  fRe    (kNSD,kNSD),
  fUe    (kNSD),
  fSpecD (kNSD),

  // 2nd order identity tensor
  fISym  (kNSD),
  fIMatx (kNSD,kNSD),

  // work spaces: 3x3 (sym & unsym) and 6x6 matrices
  fSymMatx1 (kNSD),
  fSymMatx2 (kNSD),
  fSymMatx3 (kNSD),
  fMatxCe   (kNSD,kNSD),
  fMatxSb   (kNSD,kNSD),
  fMatx1    (kNSD,kNSD),
  fMatx2    (kNSD,kNSD),
  fMatx3    (kNSD,kNSD),
  fRank4    (dSymMatrixT::NumValues(kNSD)),

  // work spaces: 9x1 arrays, 
  fFpArray (kNSD*kNSD),
  fArray1  (kNSD*kNSD),
  fArray2  (kNSD*kNSD),

  // work spaces: 9x9 matrices
  fRankIV_1 (kNSD*kNSD, kNSD*kNSD),
  fRankIV_2 (kNSD*kNSD, kNSD*kNSD),
  fRankIV_3 (kNSD*kNSD, kNSD*kNSD),
  fLHS      (kNSD*kNSD, kNSD*kNSD),

  // work spaces: arrays of 3x3 matrices
  fA           (fNumSlip),
  fB           (fNumSlip),
  fArrayOfMatx (fNumSlip),

  // average stress
  fAvgStress (kNSD)
{
  // allocate additional space for arrays of matrices
  for (int i = 0; i < fNumSlip; i++)
    {
      fA[i].Dimension(kNSD);
      fB[i].Dimension(kNSD);
      fArrayOfMatx[i].Dimension(kNSD,kNSD);
    }

  // set 2nd order unit tensors
  fIMatx.Identity();
  fISym.Identity();
}

LocalCrystalPlastFp::~LocalCrystalPlastFp() {} 

int LocalCrystalPlastFp::NumVariablesPerElement()
{
  int d_size = 0;
  int dim = dSymMatrixT::NumValues(kNSD);

  // variables per crystal per ip
  d_size += kNSD*kNSD;                       // fRotMat (const)
  d_size += kNSD*kNSD;                       // fFp_n   (t_n)
  d_size += kNSD*kNSD;                       // fFe_n   (t_n)
  d_size += kNSD*kNSD;                       // fFp     (t)
  d_size += kNSD*kNSD;                       // fFe     (t)
  d_size += fNumSlip;                        // fDGamma (t)
  d_size += dim;                             // fs_ij   (t)
  d_size += fHardening->NumberOfVariables(); // Hard Vars at t_n and t

  // total # crystal variables per element (at all IP's)
  d_size *= NumIP() * fNumGrain;

  // averaged (aggregate) stress and moduli (at all IP's)
  d_size += NumIP() * dim;                   // fsavg_ij   (t)
  d_size += NumIP() * dim * dim;             // fcavg_ijkl (t)

  return d_size;
}

int LocalCrystalPlastFp::NumberOfUnknowns() const { return kNSD*kNSD; }

const dSymMatrixT& LocalCrystalPlastFp::s_ij()
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
      fsavg_ij   = 0.0;
      //fcavg_ijkl = 0.0;    // **do not use for LS**

      // total deformation gradients
      Compute_Ftot_last_3D(fFtot_n);
      Compute_Ftot_3D(fFtot);

      for (int igrn = 0; igrn < fNumGrain; igrn++)
	{
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

	  // compute crystal Cauchy stress (elastic predictor at first iteration)
          /* use these two lines when MPS */
          //if (fFSMatSupport.StepNumber() == 0 &&        
	  //        fFSMatSupport.IterationNumber() <= -1)
          /* use these two lines when FEM */
          if (fFSMatSupport->StepNumber() >= 0 &&
              fFSMatSupport->IterationNumber() <= -1)
	     {
	       // defomation gradient
               fMatx1.SetToCombination(1., fFtot, -1., fFtot_n);
	       fFt.SetToCombination(1., fFtot_n, 100.0 / 100.0, fMatx1);

	       // elastic tensors
               fFpi.Inverse(fFp_n);
	       fFe.MultAB(fFt, fFpi);
	       fCeBar.MultATA(fFe);

	       // 2nd Piola Kirchhoff stress
               fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
	       fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

	       // Cauchy stress
               fs_ij.MultQBQT(fFe, fSBar);
               fs_ij /= fFe.Det();
	     }
	  else
             {
	       // Schmidt tensor in Bbar configuration (sample axes)
	       for (int i = 0; i < fNumSlip; i++)
                   fZ[i].MultQBQT(fRotMat, fZc[i]);

	       // compute crystal state
	       SolveCrystalState();
	  
     	       // compute crystal Cauchy stress
	       CrystalS_ij();
	  
	       // compute crystal moduli (consistent tangent)
	       //CrystalC_ijkl();  // **do not use for LS**
	     }

	  // add stress and moduli to corresponding averaged quantities
	  fsavg_ij.AddScaled(1./fNumGrain, fs_ij);
	  //fcavg_ijkl.AddScaled(1./fNumGrain, fc_ijkl); // **do not use for LS**
	}
    }

  // return averaged stress
  return fsavg_ij;
}

/*const dMatrixT& LocalCrystalPlastFp::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // return averaged moduli
  return fcavg_ijkl;
}*/

/*  **use when LS is on ** */
const dMatrixT& LocalCrystalPlastFp::c_ijkl()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // initialize averaged moduli at IP
  fcavg_ijkl = 0.0;

  // compute averaged moduli 
  for (int igrn = 0; igrn < fNumGrain; igrn++)
     {
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

	// compute consistent tangent (elastic predictor at fisrt iteration)
        if (fFSMatSupport->StepNumber() >= 1 &&
	        fFSMatSupport->IterationNumber() <= 0)
	    {
               // elastic crystal stiffness
               FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
	    }
	else
            {
	       // Schmidt tensor in Bbar configuration (sample axes)
	       for (int i = 0; i < fNumSlip; i++)
                   fZ[i].MultQBQT(fRotMat, fZc[i]);

	       // elasto-plastic crystal moduli (consistent tangent)
	       CrystalC_ijkl();  
	    }

	// add moduli to corresponding averaged quantities
	fcavg_ijkl.AddScaled(1./fNumGrain, fc_ijkl);
     }

  // return averaged moduli
  return fcavg_ijkl;
}

/* form residual for local Newton iteration */
void LocalCrystalPlastFp::FormRHS(const dArrayT& fparray, dArrayT& rhs)
{
  // recover tensor form of Fp
  Rank2FromArray9x1(fFp, fparray);

  // inverse of Fp
  fFpi.Inverse(fFp);

  // elastic right Cauchy-Green tensor 
  fCeBar.MultQTBQ(fFpi, fC);

  // Resolve Shear Stress on Slip Systems
  ResolveShearStress();

  // slip shearing rate (note: used by slip hardening classes)
  for (int i = 0; i < fNumSlip; i++)
     fDGamma[i] = fdt * fKinetics->Phi(fTau[i], i);

  // compute residual : SUM(GamDot*Z) - 1/dt*(I-Fp_n*Fp^(-1)) + pen*(detFp-1)*I
  // ... term: SUM(GamDot*Z)
  fMatx1 = 0.;
  for (int i = 0; i < fNumSlip; i++)
     fMatx1.AddScaled(fDGamma[i]/fdt, fZ[i]);

  // ... SUM(GamDot*Z) - 1/dt*(I-Fp_n*Fp^(-1))
  fMatx2.MultAB(fFp_n, fFpi);
  fMatx3.SetToCombination(1.0, fIMatx, -1.0, fMatx2);
  fMatx1.AddScaled(-1.0/fdt, fMatx3);
 
  // ... SUM(GamDot*Z) - 1/dt*(I-Fp_n*Fp^(-1)) + pen*(detFp-1)*I
  fMatx1.AddScaled(fPenalty*(fFp.Det()-1.0), fIMatx);

  // 9x1 array form of residual
  Rank2ToArray9x1(fMatx1, rhs);
}

/* form Jacobian for local Newton iteration */
void LocalCrystalPlastFp::FormLHS(const dArrayT& fparray, dMatrixT& lhs)
{
  // recover tensor form of Fp
  Rank2FromArray9x1(fFp, fparray);

  // inverse of Fp
  fFpi.Inverse(fFp);

  // elastic right Cauchy-Green tensor 
  fCeBar.MultQTBQ(fFpi, fC);

  // Resolve Shear Stress on Slip Systems
  ResolveShearStress();

  // tensor:  Aa = dTau/dCe = (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
  for (int i = 0; i < fNumSlip; i++)
     {
       fMatx1.MultAB(fMatxCe, fZ[i]);   // (CeBar*Z)
       fMatx2.MultAB(fZ[i], fMatxSb);   // (Z*SBar)
       fSymMatx1.Symmetrize(fMatx1);    // (CeBar*Z)_s
       fSymMatx2.Symmetrize(fMatx2);    // (Z*SBar)_s

       fSymMatx3.A_ijkl_B_kl(fcBar_ijkl, fSymMatx1);  // (Cijkl*(CeBar*Z)_s)
       fSymMatx2.AddScaled(0.5, fSymMatx3);           // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
       fSymMatx2.ToMatrix(fArrayOfMatx[i]);
     }

  // local Jacobian: lhs = ...
  // ... + (-2) * SUM( dGamDot/dTau * Z (x) (CeBar*Aa*Fp^(-T)) ) ...
  lhs = 0.0;
  for (int i = 0; i < fNumSlip; i++)
    {
       fMatx1.MultAB(fMatxCe, fArrayOfMatx[i]);    // (CeBar*Aa[i])
       fMatx2.MultABT(fMatx1, fFpi);               // (CeBar*Aa[i])*Fpi^T
       OuterProduct9x9(fZ[i], fMatx2, fRankIV_1);  // Z (x) (CeBar*Aa)*Fpi^T
       lhs.AddScaled(-2.0*fKinetics->DPhiDTau(fTau[i],i), fRankIV_1);
    }

  // ... + (-1/dt) * ( (Fp_n*Fp^(-1)) _(x)_ Fp^(-T) ) ...
  fMatx1.MultAB(fFp_n, fFpi);
  Rank2ToMatrixAL9x9(fMatx1, fRankIV_2);
  Rank2ToMatrixAR9x9(fFpi, fRankIV_3);
  fRankIV_1.MultAB(fRankIV_2, fRankIV_3);
  lhs.AddScaled(-1.0/fdt, fRankIV_1);

  // ... + pen*det(Fp) * ( I (x) Fp^(-T) )
  fMatx1.Transpose(fFpi);
  OuterProduct9x9(fIMatx, fMatx1, fRankIV_1);
  lhs.AddScaled(fPenalty*fFp.Det(), fRankIV_1);

  // gradient term contribution (when local, it does nothing)
  AddGradTermToLHS(lhs, fFpi);
}

void LocalCrystalPlastFp::UpdateHistory()
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
	  fFp_n = fFp;
	  fFe_n = fFe;
          fHardening->UpdateHistory();
	}
    }
}

void LocalCrystalPlastFp::ResetHistory()
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
	  fFp = fFp_n;
	  fFe = fFe_n;
          fHardening->ResetHistory();
	}
    }
}

int LocalCrystalPlastFp::NumOutputVariables() const {return kNumOutput;}

void LocalCrystalPlastFp::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void LocalCrystalPlastFp::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int group = ContinuumElement().ElementGroupNumber();
  int elem  = CurrElementNumber();
  int intpt = CurrIP();

  // load aggregate data
  LoadAggregateData(element, intpt);

  // Von Mises stress of the aggregate
  output[0] = sqrt(fSymMatx1.Deviatoric(fsavg_ij).ScalarProduct())/sqrt23;
  // cerr << " S_eq = " << output[0] << endl;
  //output[0] /= (fHardening->MaterialProperties())[1];

  // compute averaged equivalent stress
  if (elem == 0 && intpt == 0) fAvgStress = 0.0;
  fAvgStress += fsavg_ij;

  // cout << " elem = " << elem << "   intpt = " << intpt << endl;
  // cout << "    fsavg_ij = " << endl << fsavg_ij << endl;
  // cout << "    fAvgStress = " << endl << fAvgStress << endl;
  if (elem == (NumElements()-1) && intpt == (NumIP()-1))
	{
	const CommunicatorT* comm = fFSMatSupport->GroupCommunicator();
	dArrayT stress_sum(fAvgStress.Length());
	comm->Sum(fAvgStress, stress_sum);

	int total_denominator = comm->Sum(NumIP()*NumElements());
	fAvgStress.SetToScaled(1.0/total_denominator, stress_sum);
	
     cerr << " step # " << fFSMatSupport->StepNumber() 
          << "    group # " << group
          << "    S_eq_avg = " 
          << sqrt(fSymMatx1.Deviatoric(fAvgStress).ScalarProduct())/sqrt23 << endl; 	}

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
          fSpecD.PolarDecomp(fFe, fRe, fUe, false);

	  // texture: compute new crystal orientation matrix
	  fMatx1.MultAB(fRe, fRotMat);

	  // texture: compute Euler angles (in radians)
	  dArrayT& angles = fangles[igrn];
	  fLatticeOrient->RotMatrixToAngles(fMatx1, angles);
	}

      // write texture at IP/ELE
      fLatticeOrient->WriteTexture(group, elem, intpt, fNumGrain, step, fangles);
    }
}

GlobalT::SystemTypeT LocalCrystalPlastFp::TangentType() const
{
  return GlobalT::kNonSymmetric;
}

/* PROTECTED MEMBER FUNCTIONS */

void LocalCrystalPlastFp::SetSlipKinetics()
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
      throwRunTimeError("LocalCrystalPlastFp::SetSlipKinetics: Bad fKinEqnCode");
    }
  if (!fKinetics) throwMemoryError("LocalCrystalPlastFp::SetSlipKinetics");
}

void LocalCrystalPlastFp::SetSlipHardening()
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
      throwRunTimeError("LocalCrystalPlastFp::SetSlipHardening: Not implemented");
      break;

    case SlipHardening::kHard_L3:           // Haasen's model (iso)
      fHardening = new HaasenHardening(*this);
      break;

    default:
      throwRunTimeError("LocalCrystalPlastFp::SetSlipHardening: Bad fHardCode");
    }
  if (!fHardening) throwMemoryError("LocalCrystalPlastFp::SetSlipHardening");
}

void LocalCrystalPlastFp::InitializeCrystalVariables(ElementCardT& element)
{
	int elem = CurrElementNumber();

	// ... at each integration point and ...
	for (int intpt = 0; intpt < NumIP(); intpt++)
	{
	  // load aggregate data at integration point
	  LoadAggregateData(element, intpt);

	  // initialilize average stress and moduli
	  fsavg_ij = 0.0;
	  fcavg_ijkl = 0.0;

	  // ... at each crystal
	  for (int igrn = 0; igrn < fNumGrain; igrn++)
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

	      // shear rates on slip systems
	      fDGamma = 0.;

	      // crystal Cauchy stress
	      fs_ij = 0.;

	      // hardening variables
	      fHardening->InitializeHardVariables(); 
	    }
	}
}

void LocalCrystalPlastFp::LoadCrystalData(ElementCardT& element,
					int intpt, int igrain)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int stressdim   = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 5*kNSD*kNSD + fNumSlip + stressdim + fHardening->NumberOfVariables();
  int blockPerAgg = stressdim + stressdim*stressdim;
  int dex = intpt*(fNumGrain*blockPerGrn + blockPerAgg) + igrain*blockPerGrn;

  fRotMat.Set    (kNSD,kNSD, &d_array[dex             ]);   
  fFp_n.Set      (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fFe_n.Set      (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);     
  fFp.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);
  fFe.Set        (kNSD,kNSD, &d_array[dex += kNSD*kNSD]);
  fDGamma.Set    (fNumSlip,  &d_array[dex += kNSD*kNSD]);
  fs_ij.Set      (kNSD,      &d_array[dex += fNumSlip ]);     
  fHardening->LoadHardData(stressdim, dex, d_array);
}

void LocalCrystalPlastFp::LoadAggregateData(ElementCardT& element, int intpt)
{
  // fetch internal variable array
  dArrayT& d_array = element.DoubleData();
  
  // decode
  int stressdim = dSymMatrixT::NumValues(kNSD);
  int blockPerGrn = 5*kNSD*kNSD + fNumSlip + stressdim + fHardening->NumberOfVariables();
  int blockPerAgg = stressdim + stressdim*stressdim;
  int dex = intpt*(fNumGrain*blockPerGrn + blockPerAgg) + fNumGrain*blockPerGrn;

  fsavg_ij.Set   (kNSD,                &d_array[dex             ]);   
  fcavg_ijkl.Set (stressdim,stressdim, &d_array[dex += stressdim]);     
}

void LocalCrystalPlastFp::IterateOnCrystalState(bool& stateConverged, int subIncr)
{
  // initialze flag to track convergence
  stateConverged = false;

  // right Cauchy-Green tensor
  fC.MultATA(fFt);

  // initial guess for first sub-increment (uses fFe_n)
  if (subIncr == 1)
    {
      // initial ("heuristic") estimate for Fp
      InitialEstimateForFp();

      // explicit estimate for hardness
      InitialEstimateForHardening();
    }
  
  // iterate for crystal state
  int iterState = 0;
  int ierr = 0;

  while (++iterState <= fMaxIterState && !stateConverged)
     {
       try
          {
            // iter level 1: solve for Fp; Hardness = constant
            SolveForPlasticDefGradient(ierr);
            if (ierr != 0)
               {
                 if (XTAL_MESSAGES)
                    cout << " ...... failed at subIncr # " << subIncr 
			 << ";  elem # " << CurrElementNumber() 
                         << ";  IP # " << CurrIP() << endl;
	         return;
               }
      
            // iter level 2: solve for Hardness; fFp = constant
            SolveForHardening();
      
            // check convergence of state
            stateConverged = (Converged(fTolerState) && fHardening->Converged(fTolerState));
//              stateConverged = true;
         }
	  
       catch(ExceptionT::CodeT code)
	  {
            if (XTAL_MESSAGES)
               cout << " ...... failed at subIncr # " << subIncr 
	            << ";  elem # " << CurrElementNumber() 
                    << ";  IP # " << CurrIP() << endl;
	    break;
          }
     }

  // check if did not converge in max iterations
  if (!stateConverged && iterState > fMaxIterState) {
     if (XTAL_MESSAGES) {
	writeWarning("... in LocalCrystalPlastFp::IterateOnCrystalState: iters > maxIters");
        cout << " ...... failed at subIncr # " << subIncr 
	     << ";  elem # " << CurrElementNumber() 
             << ";  IP # " << CurrIP() << endl;
     }
    return;
  }
  
  // update iteration counter for state
  if (stateConverged) fIterState = max(fIterState, --iterState);
}

void LocalCrystalPlastFp::InitialEstimateForFp()
{
  // initial estimate for Fp
  fMatx1.Inverse(fFe_n);
  fFp.MultAB(fMatx1, fFt);
  fFp /= pow(fFp.Det(), 1./3.);

  // norm of Fp (to check convergence of state iterations)
  fFpNorm0 = sqrt(fFp.ScalarProduct());
}

void LocalCrystalPlastFp::InitialEstimateForHardening()
{
  // estimate for slip shearing rate (needed by slip hardening classes)
  fCeBar.MultATA(fFe_n);
  ResolveShearStress();
  for (int i = 0; i < fNumSlip; i++)
//     fDGamma[i] = fdt * fKinetics->Phi(fTau[i], i);
     fDGamma[i] = 0.0;

  // explicit update for hardening variables
  fHardening->ExplicitUpdateHard();
}

void LocalCrystalPlastFp::SolveForPlasticDefGradient(int& ierr)
{
  // 9x1 array form of Fp
  Rank2ToArray9x1(fFp, fFpArray);

  // uses "kind of" continuation method based on rate sensitivity exponent
  fKinetics->SetUpRateSensitivity();

  do {
       // current value for rate sensitivity exponent
       fKinetics->ComputeRateSensitivity();
 
       // solve for Fp
       try { fSolver->Solve(fSolverPtr, fFpArray, ierr); }
       catch(ExceptionT::CodeT code) 
           {
             if (XTAL_MESSAGES) 
		writeWarning("... in LocalCrystalPlastFp::SolveForPlasticDefGradient: exception caugth");
             fKinetics->RestoreRateSensitivity();
             throw;
           }

       // return if problems in NCLSolver
       if (ierr != 0) {
             if (XTAL_MESSAGES) 
                writeWarning("... in LocalCrystalPlastFp::SolveForPlasticDefGradient: ierr!=0");
             fKinetics->RestoreRateSensitivity();
             return;
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

void LocalCrystalPlastFp::SolveForHardening()
{
  // slip shearing rate (needed by slip hardening classes)
  fFpi.Inverse(fFp);
  fCeBar.MultQTBQ(fFpi, fC);
  ResolveShearStress();
  for (int i = 0; i < fNumSlip; i++)
     fDGamma[i] = fdt * fKinetics->Phi(fTau[i], i);

  // implicit solution for hardening variables
  //fHardening->ImplicitSolveHard();
  fHardening->ImplicitUpdateHard();
}

void LocalCrystalPlastFp::RestoreSavedSolution()
{
  fFp = fFp_save;
  fHardening->RestoreSavedSolution();
}

void LocalCrystalPlastFp::SaveCurrentSolution()
{
  fFp_save = fFp;
  fHardening->SaveCurrentSolution();
}

bool LocalCrystalPlastFp::Converged(double toler)
{
  // check convergence on Fp
  bool test = ( fabs(fFpNorm-fFpNorm0) < toler*fFpNorm0 );

  // if did not converge, reset norm
  if (!test) fFpNorm0 = fFpNorm;

  return test;
}

void LocalCrystalPlastFp::CrystalS_ij()
{
  // inverse of Fp
  fFpi.Inverse(fFp);

  // Elastic Deformation Gradient
  fFe.MultAB(fFtot, fFpi);

  // elastic rigth Cauchy-Green tensor
  fCeBar.MultATA(fFe);

  // 2nd Piola Kirchhoff Stress
  fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
  fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

  // Cauchy Stress
  fs_ij.MultQBQT(fFe, fSBar);
  fs_ij /= fFe.Det();
}

void LocalCrystalPlastFp::CrystalC_ijkl()
{
  // inverse of Fp (**use when LS is on**)
  fFpi.Inverse(fFp);

  // elastic right Cauchy-Green tensor (**use when LS is on**)
  fCeBar.MultATA(fFe);

  // Resolve Shear Stress on Slip Systems (**use when LS is on**)
  ResolveShearStress();

  // Resolve Shear Stress (fCeBar & fSBar from CrystalS_ij)
  /*fCeBar.ToMatrix(fMatxCe);
  fSBar.ToMatrix(fMatxSb);
  fMatx1.MultAB(fMatxCe, fMatxSb); 
  for (int i = 0; i < fNumSlip; i++)
     fTau[i] = dArrayT::Dot(fMatx1, fZ[i]); */

  // 2nd TERM : B(=Aa) = dTau/dCe = (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
  for (int i = 0; i < fNumSlip; i++)
    {
      fMatx1.MultAB(fMatxCe, fZ[i]);   // (CeBar*Z)
      fMatx2.MultAB(fZ[i], fMatxSb);   // (Z*SBar)
      fSymMatx1.Symmetrize(fMatx1);    // (CeBar*Z)_s
      fSymMatx2.Symmetrize(fMatx2);    // (Z*SBar)_s

      fSymMatx3.A_ijkl_B_kl(fcBar_ijkl, fSymMatx1); // (Cijkl*(CeBar*Z)_s)
      fSymMatx2.AddScaled(0.5, fSymMatx3);   // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
      fB[i].SetToScaled(1., fSymMatx2);
    }

  // compute modified Schmidt tensor : X = fLHS^(-1) * Z i  ->  fArrayOfMatx
  // where : fLHS = -2*SUM(dGamDot/dTau * Z (x) (CeBar*Aa*I)) ...
  //            ... -1/dt*(Fp_n*Fp^(-1)) _(x)_ I + pen*detFp*(I (x) I)
  ComputeSchmidtXTensor();

  // 1st TERM : A = dTau/dCe = (X*SBar)_s + 0.5*(Cijkl*(CeBar*X)_s)
  for (int i = 0; i < fNumSlip; i++)
    {
      // modified dTau/dCe
      dMatrixT& Xa = fArrayOfMatx[i];
      fMatx1.MultAB(fMatxCe, Xa);       // (CeBar*X)
      fMatx2.MultAB(Xa, fMatxSb);       // (X*SBar)
      fSymMatx1.Symmetrize(fMatx1);     // (CeBar*X)_s
      fSymMatx2.Symmetrize(fMatx2);     // (X*SBar)_s

      fSymMatx3.A_ijkl_B_kl(fcBar_ijkl, fSymMatx1); // (Cijkl*(CeBar*X)_s)
      fSymMatx2.AddScaled(0.5, fSymMatx3);   // (X*SBar)_s + 0.5*(Cijkl*(CeBar*X)_s)
      fA[i].SetToScaled(1.0, fSymMatx2);
    }

  // consistent tangent in Bbar:  CepBar = CeBar + 4*SUM(dGamDot/dTau*(A (x) B))
  for (int i = 0; i < fNumSlip; i++)
    {
      fRank4.Outer(fA[i], fB[i]);
      fcBar_ijkl.AddScaled(4.0*fKinetics->DPhiDTau(fTau[i],i), fRank4);
    }

  // gradient term contribution to moduli (when local, it does nothing)
  AddGradTermToC_ijkl();

  // Cep in current config: c_ijkl = Fe_iI*Fe_jJ*Fe_kK*Fe_lL*CepBar_IJKL
  FFFFC_3D(fc_ijkl, fcBar_ijkl, fFe);
}

void LocalCrystalPlastFp::ComputeSchmidtXTensor()
{
  // modified local Jacobian: fLHS = ...
  // ... + (-2) * SUM( dGamDot/dTau * Z (x) (CeBar*Aa*I ) ...
  fLHS = 0.0;
  for (int i = 0; i < fNumSlip; i++)
    {
       fB[i].ToMatrix(fMatx3);
       fMatx1.MultAB(fMatxCe, fMatx3);             // (CeBar*Aa[i])
       fMatx2.MultABT(fMatx1, fIMatx);             // (CeBar*Aa[i])*I
       OuterProduct9x9(fZ[i], fMatx2, fRankIV_1);  // Z (x) (CeBar*Aa)*I
       fLHS.AddScaled(-2.0*fKinetics->DPhiDTau(fTau[i],i), fRankIV_1);
    }

  // ... + (-1/dt) * ( Fp_n * Fp^(-1)) _(x)_ I ) ...
  fMatx1.MultAB(fFp_n, fFpi);
  Rank2ToMatrixAL9x9(fMatx1, fRankIV_2);
  Rank2ToMatrixAR9x9(fIMatx, fRankIV_3);
  fRankIV_1.MultAB(fRankIV_2, fRankIV_3);
  fLHS.AddScaled(-1.0/fdt, fRankIV_1);

  // ... + pen * det(Fp) * ( I (x) I )
  //fMatx1.Transpose(fIMatx);
  OuterProduct9x9(fIMatx, fIMatx, fRankIV_1);
  fLHS.AddScaled(fPenalty*fFp.Det(), fRankIV_1);

  // gradient term contribution (when local, it does nothing)
  AddGradTermToLHS(fLHS, fIMatx);

  // invert modified local Jacobian
  fRankIV_1 = MatrixInversion(fLHS);

  // compute tensor X: X[i] = fLHS^(-1) * Z[i]
  for (int i = 0; i < fNumSlip; i++)
    {
       Rank2ToArray9x1(fZ[i], fArray1);
       fRankIV_1.Multx(fArray1, fArray2);
       Rank2FromArray9x1(fArrayOfMatx[i], fArray2);
    }
}

void LocalCrystalPlastFp::ResolveShearStress()
{
  // 2nd Piola Kirchhoff Stress
  fEeBar.SetToCombination(0.5, fCeBar, -0.5, fISym);
  fSBar.A_ijkl_B_kl(fcBar_ijkl, fEeBar);

  // Resolve Shear Stress on Slip Systems
  fCeBar.ToMatrix(fMatxCe);
  fSBar.ToMatrix(fMatxSb);
  fMatx1.MultAB(fMatxCe, fMatxSb); 
  for (int i = 0; i < fNumSlip; i++)
     fTau[i] = dArrayT::Dot(fMatx1, fZ[i]);
}

void LocalCrystalPlastFp::Rank2ToArray9x1(const dMatrixT& matrix, dArrayT& array)
{
  // A -> { A } = { A00, A10, A20, A01, A11, A21, A02, A12, A22 }
  int ij = 0;
  for (int j = 0; j < kNSD; j++)
    for (int i = 0; i < kNSD; i++)
      array[ij++] = matrix(i, j);
}

void LocalCrystalPlastFp::Rank2FromArray9x1(dMatrixT& matrix, const dArrayT& array)
{
  // { A } ->  A
  int ij = 0;
  for (int j = 0; j < kNSD; j++)
    for (int i = 0; i < kNSD; i++)
      matrix(i, j) = array[ij++];
}

void LocalCrystalPlastFp::OuterProduct9x9(const dMatrixT& matrix1, 
    const dMatrixT& matrix2, dMatrixT& outer)
{
  // matrix to array: A -> { A }
  Rank2ToArray9x1(matrix1, fArray1);
  Rank2ToArray9x1(matrix2, fArray2);

  // outer product : A (x) B = { A } { B }^T
  for (int i = 0; i < kNSD*kNSD; i++)
     for (int j = 0; j < kNSD*kNSD; j++)
         outer(i, j) = fArray1[i] * fArray2[j];
}

void LocalCrystalPlastFp::Rank2ToMatrixAR9x9(const dMatrixT& matrix, dMatrixT& AR)
{
  // contructs 9x9 matrix AR : X * A = [ AR ] * { X }
  AR = 0.0;

  AR(0,0)=matrix(0,0); AR(0,3)=matrix(1,0); AR(0,6)=matrix(2,0);
  AR(1,1)=matrix(0,0); AR(1,4)=matrix(1,0); AR(1,7)=matrix(2,0);
  AR(2,2)=matrix(0,0); AR(2,5)=matrix(1,0); AR(2,8)=matrix(2,0);

  AR(3,0)=matrix(0,1); AR(3,3)=matrix(1,1); AR(3,6)=matrix(2,1);
  AR(4,1)=matrix(0,1); AR(4,4)=matrix(1,1); AR(4,7)=matrix(2,1);
  AR(5,2)=matrix(0,1); AR(5,5)=matrix(1,1); AR(5,8)=matrix(2,1);

  AR(6,0)=matrix(0,2); AR(6,3)=matrix(1,2); AR(6,6)=matrix(2,2);
  AR(7,1)=matrix(0,2); AR(7,4)=matrix(1,2); AR(7,7)=matrix(2,2);
  AR(8,2)=matrix(0,2); AR(8,5)=matrix(1,2); AR(8,8)=matrix(2,2);
}

void LocalCrystalPlastFp::Rank2ToMatrixAL9x9(const dMatrixT& matrix, dMatrixT& AL)
{
  // constructs 9x9 matrix AL : A * X = [ AL ] * { X }
  AL = 0.0;

  AL(0,0)=matrix(0,0); AL(0,1)=matrix(0,1); AL(0,2)=matrix(0,2);
  AL(1,0)=matrix(1,0); AL(1,1)=matrix(1,1); AL(1,2)=matrix(1,2);
  AL(2,0)=matrix(2,0); AL(2,1)=matrix(2,1); AL(2,2)=matrix(2,2);

  AL(3,3)=matrix(0,0); AL(3,4)=matrix(0,1); AL(3,5)=matrix(0,2);
  AL(4,3)=matrix(1,0); AL(4,4)=matrix(1,1); AL(4,5)=matrix(1,2);
  AL(5,3)=matrix(2,0); AL(5,4)=matrix(2,1); AL(5,5)=matrix(2,2);

  AL(6,6)=matrix(0,0); AL(6,7)=matrix(0,1); AL(6,8)=matrix(0,2);
  AL(7,6)=matrix(1,0); AL(7,7)=matrix(1,1); AL(7,8)=matrix(1,2);
  AL(8,6)=matrix(2,0); AL(8,7)=matrix(2,1); AL(8,8)=matrix(2,2);
}

void LocalCrystalPlastFp::AddGradTermToLHS(dMatrixT& lhs, const dMatrixT& matx)
{
  // do nothing
#pragma unused(lhs)
#pragma unused(matx)
}

void LocalCrystalPlastFp::AddGradTermToC_ijkl()
{
  // do nothing
}

/*void LocalCrystalPlastFp::InitialEstimateForHardening()
{
  // reference to hardening/kinetics material properties
  const dArrayT& propH = fHardening->MaterialProperties();
  const dArrayT& propKE = fKinetics->MaterialProperties();
  double m = propKE[0];
  double time = fFSMatSupport->Time();

  // some local tensors
  dSymMatrixT fCeDot (kNSD);

  // elasticity tensor at t_n
  fCeBar.MultATA(fFe_n);
 
  // rate : 0.5*L_v^p(CeBar)
  fMatx1.Inverse(fFtot_n);
  fMatx2.MultAB(fFt, fMatx1);                // frel
  fMatx3.Inverse(fMatx2);                    // (frel)^(-1)
  fSymMatx1.MultATA(fMatx2);                 
  fSymMatx1.PlusIdentity(-1.);
  fSymMatx1 *= 0.5 / fdt;                    // Edot=0.5*(frel^T frel - I)/dt
  fSymMatx2.MultQTBQ(fMatx3, fSymMatx1);     // Edot=frel^(-T)*Edot*frel^(-1)
  fCeDot.MultQTBQ(fFe_n, fSymMatx2);         // Cedot = 0.5*L_v^p(CeBar) = Fe^T*Edot*Fe
  //fCeDot.MultQTBQ(fFe_n, fSymMatx1);         // Cedot = 0.5*L_v^p(CeBar) = Fe^T*Edot*Fe

  // Resolve Shear Stress on Slip Systems
  ResolveShearStress();

  for (int i = 0; i < fNumSlip; i++)
    {
      fMatx1.MultAB(fMatxCe, fZ[i]);   // (CeBar*Z)
      fMatx2.MultAB(fZ[i], fMatxSb);   // (Z*SBar)
      fB[i].Symmetrize(fMatx1);        // (CeBar*Z)_s
      fA[i].Symmetrize(fMatx2);        // (Z*SBar)_s
    }

  // Form LHS
  LAdMatrixT lhs(fNumSlip);
  lhs = 0.;

  fCeDot.ToMatrix(fMatx3);
  double theta = 0.5;
       dArrayT fDGamma_n(fNumSlip);
       fDGamma_n = 0.;

  for (int i = 0; i < fNumSlip; i++)
    {
      fSymMatx3.A_ijkl_B_kl(fcBar_ijkl, fB[i]); // (Cijkl*(CeBar*Z)_s)
      fA[i].AddScaled(0.5, fSymMatx3);   // (Z*SBar)_s + 0.5*(Cijkl*(CeBar*Z)_s)
      fA[i].ToMatrix(fMatx2);

      double tmp =  2. * fdt * dArrayT::Dot(fMatx2, fMatx3);
      if (time >= 0.0) 
      //if (ftime <= fdt)
        fDGamma[i] = tmp;
      else
	fDGamma[i] = fDGamma_n[i] + theta*fabs(fDGamma_n[i])/(m*fabs(fTau[i]))*tmp;

      for (int j = 0; j < fNumSlip; j++)
	{
           fSymMatx1.SetToScaled(2.0, fB[j]);    // 2*(CeBar*Z)_s
	   fSymMatx1.ToMatrix(fMatx1);

	   // (Z*SBar)_s+0.5*(Cijkl*(CeBar*Z)_s):(2*CeBar*Z)_s
	   lhs(i,j) = dArrayT::Dot(fMatx2, fMatx1);

	   if (time >= 0.0) {
	   //if (ftime <= fdt) {
	      //if (i == j) lhs(i, j) += fHardening->HardeningModulus();
	      if (i == j) lhs(i, j) += 1.;
	   }
	   else {
	      lhs(i, j) *= theta*fabs(fDGamma_n[i])/(m*fabs(fTau[i]));
	      if (i == j) lhs(i, j) += 1.;
	   }
	}
    }

  // solve for initial estimate of dgamma
  lhs.LinearSolve(fDGamma);

  // explicit update for hardening variables
  fHardening->ExplicitUpdateHard();
}*/
