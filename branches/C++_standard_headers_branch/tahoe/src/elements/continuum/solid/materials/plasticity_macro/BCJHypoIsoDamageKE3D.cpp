/* $Id: BCJHypoIsoDamageKE3D.cpp,v 1.8 2004-07-15 08:29:14 paklein Exp $ */
#include "BCJHypoIsoDamageKE3D.h"
#include "NLCSolver.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "BCJKineticEqn.h"

using namespace Tahoe;

const double sqrt32 = sqrt(3.0/2.0);

/* printing for debugging */
const bool BCJ_DMG_MESSAGES = false;
const int IPprint = -1;

/* spatial dimensions of the problem */
const int kNSD = 3;

/* number of internal scalar variables */
const int kNumInternal = 5;
const int kNumEQValues = 6;

/* number of material properties */
const int kNumMatProp = 6;

/* element output data */
const int kNumOutput = 8;
static const char* Labels[kNumOutput] = {"EQPe","EQPh","EQXie","EQXih",
	                                 "VMISES","ALPHA","KAPPA","VVF"};

BCJHypoIsoDamageKE3D::BCJHypoIsoDamageKE3D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("BCJHypoIsoDamageKE_3D"),
  BCJHypo3D(in, support),  
  fVoidGrowthModel (NULL)
{
  // re-assigning values to base class variables
  fNumInternal = kNumInternal;      // fDEQPe, fDEQPh, fALPH, fKAPP, fDAMG
  fNumEQValues = kNumEQValues;      // fEQPe_n, fEQPe, fEQPh_n, fEQPh, fEQPDot_n, fEQPDot
  
  // re-dimensioning some arrays of base class BCJHypo3D
  fInternal_n.Dimension(fNumInternal);
  fInternal.Dimension(fNumInternal);
  fInt_save.Dimension(fNumInternal);
  fEQValues.Dimension(fNumEQValues);
  fRHS.Dimension(fNumInternal);
  fLHS.Dimension(fNumInternal,fNumInternal);
  farray.Dimension(fNumInternal);

  // allocate space for derivatives of: (i) effective stresses, (ii) factor eta
  fdEQPeDot.Dimension(fNumInternal);
  fdEQPhDot.Dimension(fNumInternal);
  fdEta.Dimension(fNumInternal);
}

BCJHypoIsoDamageKE3D::~BCJHypoIsoDamageKE3D() 
{
  delete fVoidGrowthModel;
} 

void BCJHypoIsoDamageKE3D::Initialize()
{
  // inherited
  BCJHypo3D::Initialize();

  // set void growth model
  SetVoidGrowthModel();
}

const dSymMatrixT& BCJHypoIsoDamageKE3D::s_ij()
{
  // fetch current element and int point
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // compute state, stress and moduli 
  if (MaterialSupport().RunState() == GlobalT::kFormRHS)
    {
      // reset iteration counter to check NLCSolver
      if (intpt == 0) fIterCount = 0;

      //compute 3D total deformation gradient
      Compute_Ftot_3D(fFtot);
      Compute_Ftot_last_3D(fFtot_n);

      // compute state (stress and state variables)
      SolveState();
    }

  // return cauchy stress
  return fs_ij;
}

void BCJHypoIsoDamageKE3D::FormRHS(const dArrayT& array, dArrayT& rhs)
{
  // preliminary computations
  ComputeInternalQntsRHS(array);
		
  // form residuals
  // 1. from evolution equation of Dev Cauchy stress: function G1
  rhs[0] = array[kEQXie] - fEQXieTr 
             + (3.*fmu+fMatProp[1]*(1.-fEta))*fdt*fEQPeDot;
	  
  // 2. from evolution equation of Hyd Cauchy stress: function G2
  rhs[1] = array[kEQXih] - fEQXihTr + fbulk*fdt*fEQPhDot;
	  
  // 3. from evolution equation of alpha: function G3
  rhs[2] = array[kALPH] - (1.-fEta)*fXMag;

  // 4. from evolution equation of kappa: function G4
  rhs[3] = array[kKAPP] - fInternal_n[kKAPP]  
             - (fMatProp[4]-fMatProp[3]*array[kKAPP]*array[kKAPP])*fdt*fEQPeDot
             + fMatProp[5]*fdt*array[kKAPP]*array[kKAPP];
 
  // 5. from evolution equation of Isotropic Damage: function G5
  rhs[4] = array[kDAMG] - fInternal_n[kDAMG] - (1.-array[kDAMG])*fdt*fEQPhDot;
}

void BCJHypoIsoDamageKE3D::FormLHS(const dArrayT& array, dMatrixT& lhs)
{
  // preliminary computations ("ComputeInternalQntsRHS" already called)
  ComputeInternalQntsLHS(array);

  // previous operations to evaluate local Jacobian
  fsymmatx1.SetToCombination(1., fX, -sqrt32*fMatProp[1]*(1.-fEta)*fdt*fEQPeDot, fA);
  fUnitNorm.ToMatrix(fmatx1);
  fX.ToMatrix(fmatx2);
  fUnitM.ToMatrix(fmatx3);
  fsymmatx1.ToMatrix(fmatx4);

  double cnx  = dMatrixT::Dot(fmatx1, fmatx2);
  double cmn  = dMatrixT::Dot(fmatx1, fmatx3);
  double cmxa = dMatrixT::Dot(fmatx3, fmatx4);

  // Jacobian
  // 1. d(G1)/d(EQXie),d(G1)/d(EQXih), d(G1)/d(ALPH), d(G1)/(dKAPP), d(G1)/(dDAMG)
  double c = 3.*fmu + fMatProp[1]*(1.-fEta);
  lhs(0,0) = 1. - 1./sqrt32*cnx*fdEta[0] + c*fdt*fdEQPeDot[0];
  lhs(0,1) =    - 1./sqrt32*cnx*fdEta[1] + c*fdt*fdEQPeDot[1];
  lhs(0,2) =    - 1./sqrt32*cnx*fdEta[2] + c*fdt*fdEQPeDot[2];
  lhs(0,3) =    - 1./sqrt32*cnx*fdEta[3] + c*fdt*fdEQPeDot[3];
  lhs(0,4) =    - 1./sqrt32*cnx*fdEta[4] + c*fdt*fdEQPeDot[4];

  // 2. d(G2)/d(EQXie),d(G2)/d(EQXih), d(G2)/d(ALPH), d(G2)/(dKAPP), d(G2)/(dDAMG)
  lhs(1,0) =      fbulk*fdt*fdEQPhDot[0];
  lhs(1,1) = 1. + fbulk*fdt*fdEQPhDot[1];
  lhs(1,2) =      fbulk*fdt*fdEQPhDot[2];
  lhs(1,3) =      fbulk*fdt*fdEQPhDot[3];
  lhs(1,4) =      fbulk*fdt*fdEQPhDot[4];
  
  // 3. d(G3)/d(EQXie),d(G3)/d(EQXih), d(G3)/d(ALPH), d(G3)/(dKAPP), d(G3)/(dDAMG)
  c = sqrt32*(1.-fEta)*fMatProp[1]*cmn;
  lhs(2,0) =      cmxa*fdEta[0] - c*fdt*fdEQPeDot[0];
  lhs(2,1) =      cmxa*fdEta[1] - c*fdt*fdEQPeDot[1];
  lhs(2,2) = 1. + cmxa*fdEta[2] - c*fdt*fdEQPeDot[2];
  lhs(2,3) =      cmxa*fdEta[3] - c*fdt*fdEQPeDot[3];
  lhs(2,4) =      cmxa*fdEta[4] - c*fdt*fdEQPeDot[4];

  // 4. d(G4)/d(EQXie),d(G4)/d(EQXih), d(G4)/d(ALPH), d(G4)/(dKAPP), d(G4)/(dDAMG)
  c = fMatProp[4] - fMatProp[3]*array[kKAPP]*array[kKAPP];
  lhs(3,0) =   - c*fdt*fdEQPeDot[0];
  lhs(3,1) =   - c*fdt*fdEQPeDot[1];
  lhs(3,2) =   - c*fdt*fdEQPeDot[2];
  lhs(3,3) = 1.+ 2.*(fMatProp[3]*fdt*fEQPeDot + fdt*fMatProp[5])*array[kKAPP]
               - c*fdt*fdEQPeDot[3];
  lhs(3,4) =   - c*fdt*fdEQPeDot[4];
 
  // 5. d(G5)/d(EQXie),d(G5)/d(EQXih), d(G5)/d(ALPH), d(G5)/(dKAPP), d(G5)/(dDAMG)
  lhs(4,0) =                   - (1.-array[kDAMG])*fdt*fdEQPhDot[0];
  lhs(4,1) =                   - (1.-array[kDAMG])*fdt*fdEQPhDot[1];
  lhs(4,2) =                   - (1.-array[kDAMG])*fdt*fdEQPhDot[2];
  lhs(4,3) =                   - (1.-array[kDAMG])*fdt*fdEQPhDot[3];
  lhs(4,4) = 1. + fdt*fEQPhDot - (1.-array[kDAMG])*fdt*fdEQPhDot[4];
}

void BCJHypoIsoDamageKE3D::UpdateHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // update state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // update state
      fs_ij_n     = fs_ij;
      falph_ij_n  = falph_ij;
      fInternal_n = fInternal;

      // update useful scalar values
      fEQValues[kEQPeDot_n] = fEQValues[kEQPeDot];
      fEQValues[kEQPe_n]    = fEQValues[kEQPe];
      fEQValues[kEQPh_n]    = fEQValues[kEQPh];
    }
}

void BCJHypoIsoDamageKE3D::ResetHistory()
{
  // get element
  ElementCardT& element = CurrentElement();

  // reset state at each integration point
  for (int intpt = 0; intpt < NumIP(); intpt++)
    {
      // recover local data
      LoadElementData(element, intpt);

      // reset state
      fs_ij     = fs_ij_n;
      falph_ij  = falph_ij_n;
      fInternal = fInternal_n;

      // reset useful scalar values
      fEQValues[kEQPeDot] = fEQValues[kEQPeDot_n];
      fEQValues[kEQPe]    = fEQValues[kEQPe_n];
      fEQValues[kEQPh]    = fEQValues[kEQPh_n];
    }
}

int BCJHypoIsoDamageKE3D::NumOutputVariables() const {return kNumOutput;}

void BCJHypoIsoDamageKE3D::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void BCJHypoIsoDamageKE3D::ComputeOutput(dArrayT& output)
{
  // gather element/integ point information
  ElementCardT& element = CurrentElement();
  int intpt = CurrIP();

  // load element data
  LoadElementData(element, intpt);

  // equivalent deviatoric/volumetric plastic strain
  output[0] = fEQValues[kEQPe];
  output[1] = fEQValues[kEQPh];
  
  // effective deviatoric stress and pressure
  output[2] =  fInternal[kEQXie];
  output[3] = -fInternal[kEQXih];

  // mises stress
  output[4] = sqrt(fsymmatx1.Deviatoric(fs_ij).ScalarProduct())*sqrt32;
 
  // norm of alpha, kappa and void volume fraction
  output[5] = fInternal[kALPH];
  output[6] = fInternal[kKAPP];
  output[7] = fInternal[kDAMG];

  // iter counter
  //output[8] = fIterCount;

//  if (BCJ_DMG_MESSAGES && intpt == 0 && CurrElementNumber() == 0)
  if (intpt == 0 && CurrElementNumber() == 0)
     cerr << " step # " << fFSMatSupport->StepNumber()
          << " EQPe  "  << fEQValues[kEQPe] 
          << " EQXie "  << fInternal[kEQXie] 
          << " PRESS "  << -fInternal[kEQXih] 
          << " STRESS33 " << fs_ij(2,2)
          << " VMISES " << output[4]
          << " ALPHA "  << fInternal[kALPH]
	  << " KAPPA "  << fInternal[kKAPP]
          << " VVF   "  << fInternal[kDAMG]
          << " ITERS "  << fIterCount << endl;
}

/* PROTECTED MEMBER FUNCTIONS */

void BCJHypoIsoDamageKE3D::InitializeVariables(ElementCardT& element)
{
      // initialize state at each integration point
      for (int intpt = 0; intpt < NumIP(); intpt++)
	{
	  // load element data
	  LoadElementData(element, intpt);

	  // Cauchy stress
          fs_ij_n = 0.;
          fs_ij   = 0.;
	      
          // backstress
          falph_ij_n = 0.;
          falph_ij   = 0.;

	  // scalar internal variables
	  fInternal_n = 0.;
	  fInternal   = 0.;

	  // initial damage
	  fInternal_n[kDAMG] = fDamg0;
	  fInternal[kDAMG]   = fDamg0;

	  // other scalar values
	  fEQValues = 0.;

	  // tangent moduli
	  fc_ijkl = 0.;
	}
}

void BCJHypoIsoDamageKE3D::SetVoidGrowthModel()
{
  // read initial damage
  fInput >> fDamg0;
  
  // read strain rate sensitivity exponent (used in VGModel)
  fInput >> fm;

  // read void growth model code
  fInput >> fVGMCode;

  // select kinetic eqn model
  switch(fVGMCode)
    {
    case VoidGrowthModelImp::kCocks:
      fVoidGrowthModel = new CocksVGModel(*this);
      break;
      
    case VoidGrowthModelImp::kDuvaCrow:
      fVoidGrowthModel = new DuvaCrowVGModel(*this);
      break;
      
    case VoidGrowthModelImp::kSofronis:
      fVoidGrowthModel = new SofronisVGModel(*this);
      break;
      
    default:
      throwRunTimeError("BCJHypoIsoDamageKE3D::SetVoidGrowthModel: Bad fVGMCode");
    }
  if (!fVoidGrowthModel) throwMemoryError("BCJHypoIsoDamageKE3D::SetVoidGrowthModel");

  // set rate sensitivity exponent in VGModel class
  fVoidGrowthModel->SetRateSensitivity(fm);
}

void BCJHypoIsoDamageKE3D::IntegrateConstitutiveEqns(bool& converged, int subIncr,
		                                   int totSubIncrs)
{
  // initialize flag to track convergence
  converged = true;
  
  // step 1. polar decomposition of relative deformation gradient
  PolarDecomposition();

  // step 2. incremental total strain
  IncrementalStrain();

  // step 3. rotate tensorial state variables at t_n to current configuration
  RotateTensorialVariables();

  // step 4. trial stresses
  ElasticTrialStress();

  // check for inelastic process (note: uses deviatoric part)
  if ( fEQXieTr > (1.+1.e-6)*fKineticEqn->h(fEQValues[kEQPeDot_n],fInternal_n[kKAPP])
	&& fFSMatSupport->IterationNumber() > -1 )
    {
      // step 5. forward gradient estimate
      if (subIncr == 1) ForwardGradientEstimate();
      
      // step 6. solve for state variables
      Solve(converged);
      if (!converged) {
         if (BCJ_DMG_MESSAGES)
            cout << " Did not converged at elem # " << CurrElementNumber()
                 << ";  IP # " << CurrIP() << ";  subIncr/totSubIncrs = "
	         << subIncr << "/" << totSubIncrs << endl;
	 return;
      }
      if (subIncr < totSubIncrs) return;

      // step 7. compute Cauchy stress and backstress
      UpdateStresses();

      // step 8. elasto-plastic moduli
      TangentModuli();
    }
  else  // elastic process
    {
      // step 5. reset values
      fInternal = fInternal_n;
      fInternal[kEQXie] = fEQXieTr;
      fInternal[kEQXih] = fEQXihTr;
      if (subIncr < totSubIncrs) return;

      falph_ij  = falph_ij_n;
      fEQValues[kEQPeDot] = fEQValues[kEQPeDot_n];
      fEQValues[kEQPe]    = fEQValues[kEQPe_n];
      fEQValues[kEQPh]    = fEQValues[kEQPh_n];
     
      // step 6. Cauchy stress
      fs_ij = fSigTr;

      // step 7. elastic moduli
      ElasticModuli(fmu, fbulk);
    }
}

void BCJHypoIsoDamageKE3D::ElasticTrialStress()
{
  // trial stress (elastic predictor)
  fSigTr.SetToCombination(1., fsigma_n, 2.*fmu, fDE, 
			  flambda*fDE.Trace(), fISym);

  // deviatoric trial stress tensor and trial volumetric stress
  fSigTrDev.Deviatoric(fSigTr);
  fEQXihTr = fSigTr.Trace() / 3.;

  // trial equivalent (deviatoric) stress
  fXiTr.SetToCombination(1., fSigTrDev, -2./3., falpha_n);
  fEQXieTr = sqrt32*sqrt(fXiTr.ScalarProduct());
}

void BCJHypoIsoDamageKE3D::UpdateStresses()
{
  // preliminary computations
  ComputeInternalQntsRHS(fInternal);
 
  // current equivalent deviatoric/hydrostatic stresses
  //double EQXie = fEQXieTr - (3*fmu+fMatProp[1]*(1.-fEta))*fdt*fEQPeDot;
  //double EQXih = fEQXihTr - fbulk*fdt*fEQPhDot;

  // current equivalent deviatoric/hydrostatic/total plastic strains
  fEQValues[kEQPe]    = fEQValues[kEQPe_n] + fdt*fEQPeDot;
  fEQValues[kEQPh]    = fEQValues[kEQPh_n] + fdt*fEQPhDot;
  fEQValues[kEQPeDot] = fEQPeDot;

  // radial factor (deviatoric)
  fradial = fInternal[kEQXie]/fEQXieTr;

  // backstress
  falph_ij.SetToScaled((1.-fEta), fX);

  // Cauchy stress
  fs_ij.SetToCombination(fradial, fXiTr, 2./3., falph_ij, fInternal[kEQXih], fISym);
}

/* Consisten Tangent based on Updated Lagrangian Procedure */
void BCJHypoIsoDamageKE3D::TangentModuli()
{
  // preliminary computations
  ComputeInternalQntsLHS(fInternal);

  FormLHS(fInternal, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  ffactor = fradial + fMatProp[1]*(1.-fEta)*fdt*fEQPeDot/fEQXieTr;
  double gamma  = (1.-ffactor)/(1.-fradial);

  fUnitNorm.ToMatrix(fmatx1);
  fUnitM.ToMatrix(fmatx2);
  double cmn = dMatrixT::Dot(fmatx1, fmatx2);

  double tmpa1 = (2.*fmu)*sqrt32;
  double tmpa2 = (2.*fmu)*1.5*(1.-gamma)*cmn;
  double tmpc = 1.5*(2.*fmu)/fXMag*(ffactor-fradial);

  double a1 = tmpa1*Jaci(0,0) + tmpa2*Jaci(0,2);
  double a2 = tmpa1*Jaci(1,0) + tmpa2*Jaci(1,2);
  double a3 = tmpa1*Jaci(2,0) + tmpa2*Jaci(2,2);
  double a4 = tmpa1*Jaci(3,0) + tmpa2*Jaci(3,2);
  double a5 = tmpa1*Jaci(4,0) + tmpa2*Jaci(4,2);

  double b1 = fbulk*Jaci(0,1);
  double b2 = fbulk*Jaci(1,1);
  double b3 = fbulk*Jaci(2,1);
  double b4 = fbulk*Jaci(3,1);
  double b5 = fbulk*Jaci(4,1);

  double c1 = tmpc*Jaci(0,2);
  double c2 = tmpc*Jaci(1,2);
  double c3 = tmpc*Jaci(2,2);
  double c4 = tmpc*Jaci(3,2);
  double c5 = tmpc*Jaci(4,2);

  double a6 = a1*fdEta[0]+a2*fdEta[1]+a3*fdEta[2]+a4*fdEta[3]+a5*fdEta[4];
  double b6 = b1*fdEta[0]+b2*fdEta[1]+b3*fdEta[2]+b4*fdEta[3]+b5*fdEta[4];
  double c6 = c1*fdEta[0]+c2*fdEta[1]+c3*fdEta[2]+c4*fdEta[3]+c5*fdEta[4];

  // elastic-like terms: 2*mu*f*(I) +(k*g-2/3*mu*f)*(1(x)1)
  ElasticModuli(fmu*ffactor, b2);

  // inelastic part: term  beta1*(N(x)N)
  double beta = - gamma * (2.*fmu*fradial + 2./3.*a6*cmn*fXMag - a1/sqrt32);
  fRank4.Outer(fUnitNorm, fUnitNorm);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta2*(A(x)A)
  beta = - (1.-ffactor)*c6*fEQXieTr*fEQXieTr;
  fRank4.Outer(fA, fA);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta3*(I(x)N)
  beta = a2;
  fRank4.Outer(fISym, fUnitNorm);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta4*(N(x)I)
  beta = - gamma * (2./3.*cmn*fXMag*b6 - b1/sqrt32);
  fRank4.Outer(fUnitNorm, fISym);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta5*(I(x)A)
  beta = sqrt32*c2*fEQXieTr;
  fRank4.Outer(fISym, fA);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta6*(A(x)I)
  beta = - 1./sqrt32*(1.-ffactor)*b6*fEQXieTr;
  fRank4.Outer(fA, fISym);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta7*(N(x)A)
  beta = - gamma * (cmn*fXMag*c6/sqrt32 - c1) * fEQXieTr;
  fRank4.Outer(fUnitNorm, fA);
  fc_ijkl.AddScaled(beta, fRank4);

  // inelastic part: term  beta8*(A(x)N)
  beta = - 1./sqrt32*(1.-ffactor)*a6*fEQXieTr;
  fRank4.Outer(fA, fUnitNorm);
  fc_ijkl.AddScaled(beta, fRank4);
}

void BCJHypoIsoDamageKE3D::ForwardGradientEstimate()
{
  // residual at t_n (should be different from zero)
  FormRHS(fInternal_n, fRHS);

  // inverse of Jacobian at t_n (MatrixInversion changes fLHS)
  FormLHS(fInternal_n, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  // forward gradient estimate
  Jaci.Multx(fRHS, farray);
  fInternal.SetToCombination(1., fInternal_n, -1., farray);

  // checks on initial guess for deviatoric effective stress and damage
  if (fInternal[kEQXie] <= 0.0)   fInternal[kEQXie] = 1.e-10;
  if (fInternal[kALPH]  <= 0.0)   fInternal[kALPH] = 1.e-10;
  if (fInternal[kKAPP]  <= 0.0)   fInternal[kKAPP] = 1.e-10;
  if (fInternal[kDAMG]  <= 0.0)   fInternal[kDAMG] = 1.e-10;
  if (fInternal[kDAMG]  >= 0.999) fInternal[kDAMG] = 0.999;
}

bool BCJHypoIsoDamageKE3D::IsSolnVariableNegative()
{
  // test for negative values of internal variables
  // note: hydrostatic stress EQXih can be negative
  int i = 0;
  bool isNegative = false;
  while (i < fNumInternal && !isNegative)
    {
      if (i == kEQXih) break;
      isNegative = (fInternal[i] < 0.0);
      i++;
    }

  // check for zero values of internal variables
  if (!isNegative) {
     if (fInternal[kEQXie] <= 1.e-10)       fInternal[kEQXie] = 1.e-10;
     if (fabs(fInternal[kEQXih]) <= 1.e-10) fInternal[kEQXih] = 1.e-10;
     if (fInternal[kALPH] <= 1.e-10)        fInternal[kALPH] = 1.e-10;
     if (fInternal[kKAPP] <= 1.e-10)        fInternal[kKAPP] = 1.e-10;
     if (fInternal[kDAMG] <= 1.e-10)        fInternal[kDAMG] = 1.e-10;
     if (fInternal[kDAMG] >= 0.999)         fInternal[kDAMG] = 0.999;
  }

  return isNegative;
}

/* PRIVATE MEMBER FUNCTIONS */

void BCJHypoIsoDamageKE3D::ComputeInternalQntsRHS(const dArrayT& array)
{
  // effective strain rates of porous material
  // ... coefficients A1 & A2 from void growth model
  fVoidGrowthModel->ACoefficients(array[kDAMG], fA1Dmg, fA2Dmg);
  if (array[kEQXih] <= 0.0) { fA2Dmg = 0.0; }

  // ... equivalent stress of porous material
  fEQXi = sqrt(fA1Dmg*array[kEQXie]*array[kEQXie] + fA2Dmg*array[kEQXih]*array[kEQXih]);
  if (fabs(fEQXi) <= 1.e-6) fEQXi = 1.e-6;
 
  // ... equivalent plastic strain rate of porous material
  fEQPDot = fKineticEqn->f(fEQXi, array[kKAPP]);
  
  // ... deviatoric and hydrostatic effective strain rates
  fEQPeDot = fA1Dmg * array[kEQXie]/fEQXi * fEQPDot;
  fEQPhDot = fA2Dmg * array[kEQXih]/fEQXi * fEQPDot;
 
  // factor eta
  fEta = fMatProp[0]*array[kALPH]*fdt*fEQPeDot + fMatProp[2]*array[kALPH]*fdt;
  fEta = fEta/(1.+fEta);

  // "new" trial equivalente stress and unit normal in stress space
  fXiTr.SetToCombination(1., fSigTrDev, -2./3.*(1.-fEta), falpha_n);
  fEQXieTr = sqrt32*sqrt(fXiTr.ScalarProduct());
  fUnitNorm.SetToScaled(sqrt32/fEQXieTr, fXiTr);

  // factor (tensor) in integrated eqn for backstress: alpha=(1-eta)*X
  fX.SetToCombination(1., falpha_n, sqrt32*fMatProp[1]*fdt*fEQPeDot, fUnitNorm);
  fXMag = sqrt(fX.ScalarProduct());
  if (fXMag <= 1.e-16) fXMag = 1.;  // to avoid NaN in fUnitM in first incr/iter
}

void BCJHypoIsoDamageKE3D::ComputeInternalQntsLHS(const dArrayT& array)
{
  // preliminary computations
  double beta     = fKineticEqn->DfDsigma(fEQXi, array[kKAPP]);
  double coefBeta = (beta*fEQXi-fEQPDot)/fEQXi;
  double ratioE   = array[kEQXie]/fEQXi;
  double ratioH   = array[kEQXih]/fEQXi;

  // derivatives dA1/d() & dA2/d() from void growth model
  fVoidGrowthModel->ADerivCoefficients(array[kDAMG], fdA1Dmg, fdA2Dmg);
  if (array[kEQXih] <= 0.0) { fdA1Dmg = 0.0; fdA2Dmg = 0.0; }

  // derivative of (total) equivalent plastic strain increment wrt damage
  double dEQXidDAMG = 0.5*fEQXi*(ratioE*ratioE*fdA1Dmg + ratioH*ratioH*fdA2Dmg);

  // derivatives of deviatoric equivalent plastic strain rate wrt primary variables
  fdEQPeDot[0] = coefBeta*fA1Dmg*fA1Dmg*ratioE*ratioE + fA1Dmg*fEQPDot/fEQXi;
  fdEQPeDot[1] = coefBeta*fA1Dmg*fA2Dmg*ratioE*ratioH;
  fdEQPeDot[2] = 0.;
  fdEQPeDot[3] = fA1Dmg*ratioE * fKineticEqn->DfDs(fEQXi, array[kKAPP]);
  fdEQPeDot[4] = coefBeta*fA1Dmg*ratioE*dEQXidDAMG + ratioE*fEQPDot*fdA1Dmg;

  // derivatives of volumetric equivalent plastic strain rate wrt primary variables
  fdEQPhDot[0] = coefBeta*fA1Dmg*fA2Dmg*ratioE*ratioH;
  fdEQPhDot[1] = coefBeta*fA2Dmg*fA2Dmg*ratioH*ratioH + fA2Dmg*fEQPDot/fEQXi;
  fdEQPhDot[2] = 0.;
  fdEQPhDot[3] = fA2Dmg*ratioH * fKineticEqn->DfDs(fEQXi, array[kKAPP]);
  fdEQPhDot[4] = coefBeta*fA2Dmg*ratioH*dEQXidDAMG + ratioH*fEQPDot*fdA2Dmg;

  // derivatives of factor Eta wrt primary variables
  double c = (1.-fEta)*(1.-fEta)*(fMatProp[0]*array[kALPH]);
  fdEta[0] = c*fdt*fdEQPeDot[0];
  fdEta[1] = c*fdt*fdEQPeDot[1];
  fdEta[2] = (1.-fEta)*(1.-fEta)*(fMatProp[0]*fdt*fEQPeDot+fMatProp[2]*fdt);
  fdEta[3] = c*fdt*fdEQPeDot[3];
  fdEta[4] = c*fdt*fdEQPeDot[4];

  // unit normal M: M=X/||X||
  fUnitM.SetToScaled(1./fXMag, fX);

  // derivative factor (tensor) for unit normal: dN/d()=A*dEta/d()
  falpha_n.ToMatrix(fmatx1);
  fUnitNorm.ToMatrix(fmatx2);
  fsymmatx1.SetToCombination(1., falpha_n, -dArrayT::Dot(fmatx1,fmatx2), fUnitNorm);
  fA.SetToScaled(1./(sqrt32*fEQXieTr), fsymmatx1);
}
