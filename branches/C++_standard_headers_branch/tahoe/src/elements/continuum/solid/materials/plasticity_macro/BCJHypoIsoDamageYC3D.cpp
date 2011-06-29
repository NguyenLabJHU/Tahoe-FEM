/* $Id: BCJHypoIsoDamageYC3D.cpp,v 1.9 2004-07-15 08:29:14 paklein Exp $ */
#include "BCJHypoIsoDamageYC3D.h"
#include "NLCSolver.h"
#include "ElementCardT.h"

#include "Utils.h"
#include "BCJKineticEqn.h"

using namespace Tahoe;

const double sqrt32 = sqrt(3.0/2.0);

/* printing for debugging */
const bool BCJ_DMG_MESSAGES = false;
const int IPprint = -1;

/* lower bound for damage or vvf */
const double DMG_LWB_LIMIT = 1.e-6;

/* spatial dimensions of the problem */
const int kNSD = 3;

/* number of internal scalar variables */
const int kNumInternal = 5;
const int kNumEQValues = 10;

/* number of material properties */
const int kNumMatProp = 6;

/* element output data */
const int kNumOutput = 8;
static const char* Labels[kNumOutput] = {"EQPe","EQPh","EQXie","EQXih",
                                         "VMISES","ALPHA","KAPPA","VVF"};

BCJHypoIsoDamageYC3D::BCJHypoIsoDamageYC3D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("BCJHypoIsoDamageYC_3D"),
  BCJHypo3D(in, support),  
  fVoidGrowthModel (NULL)
{
  // re-assigning values to base class variables
  fNumInternal = kNumInternal;      // fDEQPe, fDEQPh, fALPH, fKAPP, fDAMG
  fNumEQValues = kNumEQValues;      // fEQPe_n, fEQPe, fEQPh_n, fEQPh, fEQP_n,
                                    // fEQP, fEQXie_n, fEQXie, fEQXih_n, fEQXih 
  
  // re-dimensioning some arrays of base class BCJHypo3D
  fInternal_n.Dimension(fNumInternal);
  fInternal.Dimension(fNumInternal);
  fInt_save.Dimension(fNumInternal);
  fEQValues.Dimension(fNumEQValues);
  fRHS.Dimension(fNumInternal);
  fLHS.Dimension(fNumInternal,fNumInternal);
  farray.Dimension(fNumInternal);

  // allocate space for derivatives of effective stresses
  fdEQXie.Dimension(fNumInternal);
  fdEQXih.Dimension(fNumInternal);
}

BCJHypoIsoDamageYC3D::~BCJHypoIsoDamageYC3D() 
{
  delete fVoidGrowthModel;
} 

void BCJHypoIsoDamageYC3D::Initialize()
{
  // inherited
  BCJHypo3D::Initialize();

  // set void growth model
  SetVoidGrowthModel();
}

const dSymMatrixT& BCJHypoIsoDamageYC3D::s_ij()
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
      if (CurrIP() == 0) fIterCount = 0;

      //compute 3D total deformation gradient
      Compute_Ftot_3D(fFtot);
      Compute_Ftot_last_3D(fFtot_n);

      // compute state (stress and state variables)
      SolveState();
    }

  // return cauchy stress
  return fs_ij;
}

void BCJHypoIsoDamageYC3D::FormRHS(const dArrayT& array, dArrayT& rhs)
{
  // preliminary computations
  ComputeInternalQntsRHS(array);
		
  // form residuals
  // 1. from evolution equation of Dev Cauchy stress: function G1
  rhs[0] = fEQValues[kEQXie] - fEQXieTr 
             + (3.*fmu + fMatProp[1]*(1.-fEta))*array[kDEQPe];
	  
  // 2. from evolution equation of Hyd Cauchy stress: function G2
  rhs[1] = fEQValues[kEQXih] - fEQXihTr + fbulk*array[kDEQPh];
	  
  // 3. from evolution equation of alpha: function G3
  rhs[2] = array[kALPH] - (1.-fEta)*fXMag;

  // 4. from evolution equation of kappa: function G4
  rhs[3] = array[kKAPP] - fInternal_n[kKAPP]  
             - (fMatProp[4]-fMatProp[3]*array[kKAPP]*array[kKAPP])*array[kDEQPe]
             + fMatProp[5]*fdt*array[kKAPP]*array[kKAPP];
 
  // 5. from evolution equation of Isotropic Damage: function G5
  rhs[4] = array[kDAMG] - fInternal_n[kDAMG] - (1.-array[kDAMG])*array[kDEQPh];
}

void BCJHypoIsoDamageYC3D::FormLHS(const dArrayT& array, dMatrixT& lhs)
{
  // preliminary computations ("ComputeInternalQntsRHS" already called)
  ComputeInternalQntsLHS(array);

  // previous operations to evaluate local Jacobian
  fsymmatx1.SetToCombination(1., fX, -sqrt32*fMatProp[1]*(1.-fEta)*array[kDEQPe], fA);
  fUnitNorm.ToMatrix(fmatx1);
  fX.ToMatrix(fmatx2);
  fUnitM.ToMatrix(fmatx3);
  fsymmatx1.ToMatrix(fmatx4);

  double cnx  = dMatrixT::Dot(fmatx1, fmatx2);
  double cmn  = dMatrixT::Dot(fmatx1, fmatx3);
  double cmxa = dMatrixT::Dot(fmatx3, fmatx4);

  // Jacobian
  // 1. d(G1)/d(DEQPe),d(G1)/d(DEQPh), d(G1)/d(ALPH), d(G1)/(dKAPP), d(G1)/(dDAMG)
  lhs(0,0) = fdEQXie[0] - 1./sqrt32*fDEtaDDEQP*cnx + (3.*fmu+fMatProp[1]*(1.-fEta));
  lhs(0,1) = fdEQXie[1];
  lhs(0,2) = fdEQXie[2] - 1./sqrt32*fDEtaDALPH*cnx;
  lhs(0,3) = fdEQXie[3];
  lhs(0,4) = fdEQXie[4];

  // 2. d(G2)/d(DEQPe),d(G2)/d(DEQPh), d(G2)/d(ALPH), d(G2)/(dKAPP), d(G2)/(dDAMG)
  lhs(1,0) = fdEQXih[0];
  lhs(1,1) = fdEQXih[1] + fbulk;
  lhs(1,2) = fdEQXih[2];
  lhs(1,3) = fdEQXih[3];
  lhs(1,4) = fdEQXih[4];
  
  // 3. d(G3)/d(DEQPe),d(G3)/d(DEQPh), d(G3)/d(ALPH), d(G3)/(dKAPP), d(G3)/(dDAMG)
  lhs(2,0) = fDEtaDDEQP*cmxa - sqrt32*(1.-fEta)*fMatProp[1]*cmn;
  lhs(2,1) = 0.0;
  lhs(2,2) = 1. + cmxa*fDEtaDALPH;
  lhs(2,3) = 0.0;
  lhs(2,4) = 0.0;

  // 4. d(G4)/d(DEQPe),d(G4)/d(DEQPh), d(G4)/d(ALPH), d(G4)/(dKAPP), d(G4)/(dDAMG)
  lhs(3,0) = -fMatProp[4] + fMatProp[3]*array[kKAPP]*array[kKAPP];
  lhs(3,1) = 0.0;
  lhs(3,2) = 0.0;
  lhs(3,3) = 1. + 2.*(fMatProp[3]*array[kDEQPe] + fMatProp[5]*fdt) * array[kKAPP];
  lhs(3,4) = 0.0;
 
  // 5. d(G5)/d(DEQPe),d(G5)/d(DEQPh), d(G5)/d(ALPH), d(G5)/(dKAPP), d(G5)/(dDAMG)
  lhs(4,0) = 0.0;
  lhs(4,1) = - (1.-array[kDAMG]);
  lhs(4,2) = 0.0;
  lhs(4,3) = 0.0;
  lhs(4,4) = 1. + array[kDEQPh];
}

void BCJHypoIsoDamageYC3D::UpdateHistory()
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
      fEQValues[kEQXie_n] = fEQValues[kEQXie];
      fEQValues[kEQXih_n] = fEQValues[kEQXih];
      fEQValues[kEQP_n]   = fEQValues[kEQP];
      fEQValues[kEQPe_n]  = fEQValues[kEQPe];
      fEQValues[kEQPh_n]  = fEQValues[kEQPh];
    }
}

void BCJHypoIsoDamageYC3D::ResetHistory()
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
      fEQValues[kEQXie] = fEQValues[kEQXie_n];
      fEQValues[kEQXih] = fEQValues[kEQXih_n];
      fEQValues[kEQP]   = fEQValues[kEQP_n];
      fEQValues[kEQPe]  = fEQValues[kEQPe_n];
      fEQValues[kEQPh]  = fEQValues[kEQPh_n];
    }
}

int BCJHypoIsoDamageYC3D::NumOutputVariables() const {return kNumOutput;}

void BCJHypoIsoDamageYC3D::OutputLabels(ArrayT<StringT>& labels) const
{
  // allocate space for labels
  labels.Dimension(kNumOutput);

  // copy labels
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void BCJHypoIsoDamageYC3D::ComputeOutput(dArrayT& output)
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
  output[2] =  fEQValues[kEQXie];
  output[3] = -fEQValues[kEQXih];

  // mises stress
  output[4] = sqrt(fsymmatx1.Deviatoric(fs_ij).ScalarProduct())*sqrt32;
  //
  // norm of alpha, kappa and void volume fraction
  output[5] = fInternal[kALPH];
  output[6] = fInternal[kKAPP];
  output[7] = fInternal[kDAMG];

  // iter counter
  //output[8] = fIterCount;

  if (BCJ_DMG_MESSAGES && intpt == 0 && CurrElementNumber() == 0)
     cerr << " step # " << fFSMatSupport->StepNumber()
          << " EQPe  "  << fEQValues[kEQPe] 
          << " EQXie "  << fEQValues[kEQXie]
          << " PRESS "  << -fEQValues[kEQXih]
          << " STRESS33 " << fs_ij(2,2)
          << " VMISES " << output[4]
          << " ALPHA "  << fInternal[kALPH]
	  << " KAPPA "  << fInternal[kKAPP]
          << " VVF   "  << fInternal[kDAMG]
	  << " ITERS "  << fIterCount << endl;
}

/* PROTECTED MEMBER FUNCTIONS */

void BCJHypoIsoDamageYC3D::InitializeVariables(ElementCardT& element)
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

void BCJHypoIsoDamageYC3D::SetVoidGrowthModel()
{
  // read initial damage
  fInput >> fDamg0;
  if (fDamg0 <= DMG_LWB_LIMIT) fDamg0 = DMG_LWB_LIMIT;

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
      throwRunTimeError("BCJHypoIsoDamageYC3D::SetVoidGrowthModel: Bad fVGMCode");
    }
  if (!fVoidGrowthModel) throwMemoryError("BCJHypoIsoDamageYC3D::SetVoidGrowthModel");

  // set rate sensitivity exponent in VGModel class
  fVoidGrowthModel->SetRateSensitivity(fm);
}

void BCJHypoIsoDamageYC3D::IntegrateConstitutiveEqns(bool& converged, int subIncr,
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

  // check for inelastic process (note: use deviatoric part for this check)
  if ( fEQXieTr > (1.+1.e-6)*fKineticEqn->h(((fabs(fdt) > kSmall) ? fInternal_n[kDEQPe]/fdt : 0.0),fInternal_n[kKAPP])
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
      if (subIncr < totSubIncrs) return;

      falph_ij  = falph_ij_n;
      fEQValues[kEQXie] = fEQXieTr;
      fEQValues[kEQXih] = fEQXihTr;
      fEQValues[kEQP]   = fEQValues[kEQP_n];
      fEQValues[kEQPe]  = fEQValues[kEQPe_n];
      fEQValues[kEQPh]  = fEQValues[kEQPh_n];
     
      // step 6. Cauchy stress
      fs_ij = fSigTr;

      // step 7. elastic moduli
      ElasticModuli(fmu, fbulk);
    }
}

void BCJHypoIsoDamageYC3D::ElasticTrialStress()
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

void BCJHypoIsoDamageYC3D::UpdateStresses()
{
  // preliminary computations
  ComputeInternalQntsRHS(fInternal);
 
  // current equivalent deviatoric/hydrostatic stresses
  // fEQValues[kEQSig] = fKineticEqn->h(fInternal[kDEQP]/fdt, fInternal[kKAPP]);
  fEQValues[kEQXie] = fEQXieTr - (3*fmu+fMatProp[1]*(1.-fEta))*fInternal[kDEQPe];
  fEQValues[kEQXih] = fEQXihTr - fbulk*fInternal[kDEQPh];

  // current equivalent deviatoric/hydrostatic/total plastic strains
  fEQValues[kEQPe] = fEQValues[kEQPe_n] + fInternal[kDEQPe];
  fEQValues[kEQPh] = fEQValues[kEQPh_n] + fInternal[kDEQPh];
  fEQValues[kEQP]  = fEQValues[kEQP_n]  + fDEQP;

  // radial factor (deviatoric)
  fradial = fEQValues[kEQXie]/fEQXieTr;

  // backstress
  falph_ij.SetToScaled((1.-fEta), fX);

  // Cauchy stress
  fs_ij.SetToCombination(fradial, fXiTr, 2./3., falph_ij, fEQValues[kEQXih], fISym);
}

/* Consisten Tangent based on Updated Lagrangian Procedure */
void BCJHypoIsoDamageYC3D::TangentModuli()
{
  // preliminary computations
  //ComputeInternalQntsRHS(fInternal);
  //ComputeInternalQntsLHS(fInternal);

  FormLHS(fInternal, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  ffactor = fradial + fMatProp[1]*(1.-fEta)*fInternal[kDEQPe]/fEQXieTr;
  
  double tmpa = sqrt32*(2.*fmu);
  double tmpc = 1.5*(2.*fmu)/fXMag*(ffactor-fradial);

  double a1 = tmpa*Jaci(0,0);
  double a2 = tmpa*Jaci(1,0);
  double a3 = tmpa*Jaci(2,0);

  double b1 = fbulk*Jaci(0,1);
  double b2 = fbulk*Jaci(1,1);
  double b3 = fbulk*Jaci(2,1);

  double c1 = tmpc*Jaci(0,2);
  double c2 = tmpc*Jaci(1,2);
  double c3 = tmpc*Jaci(2,2);

  double a6 = a1*fDEtaDDEQP + a3*fDEtaDALPH;
  double b6 = b1*fDEtaDDEQP + b3*fDEtaDALPH;
  double c6 = c1*fDEtaDDEQP + c3*fDEtaDALPH;

  // elastic-like terms: 2*mu*f*(I) +(k-2/3*mu*f)*(1(x)1)
  ElasticModuli(fmu*ffactor, fbulk*(1.-b2));

  // inelastic part: term  -beta1*(N(x)N)
  double beta = 2.*fmu*ffactor - 2.*fmu*(1.-sqrt32*a1);
  fRank4.Outer(fUnitNorm, fUnitNorm);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta2*(A(x)A)
  beta = (1.-ffactor)*c6*fEQXieTr*fEQXieTr;
  fRank4.Outer(fA, fA);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta3*(I(x)N)
  beta = fbulk*a2;
  fRank4.Outer(fISym, fUnitNorm);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta4*(N(x)I)
  beta = 2.*sqrt32*fmu*b1;
  fRank4.Outer(fUnitNorm, fISym);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta5*(I(x)A)
  beta = sqrt32*fbulk*c2*fEQXieTr;
  fRank4.Outer(fISym, fA);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta6*(A(x)I)
  beta = 1./sqrt32*(1.-ffactor)*b6*fEQXieTr;
  fRank4.Outer(fA, fISym);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta7*(N(x)A)
  beta = 3.*fmu*c1*fEQXieTr;
  fRank4.Outer(fUnitNorm, fA);
  fc_ijkl.AddScaled(-beta, fRank4);

  // inelastic part: term  -beta8*(A(x)N)
  beta = 1./sqrt32*(1.-ffactor)*a6*fEQXieTr;
  fRank4.Outer(fA, fUnitNorm);
  fc_ijkl.AddScaled(-beta, fRank4);
}

void BCJHypoIsoDamageYC3D::ForwardGradientEstimate()
{
  // residual at t_n (they should be different from zero)
  FormRHS(fInternal_n, fRHS);

  // inverse of Jacobian at t_n (MatrixInversion changes fLHS)
  FormLHS(fInternal_n, fLHS);
  dMatrixT Jaci = MatrixInversion(fLHS);

  // forward gradient estimate
  Jaci.Multx(fRHS, farray);
  fInternal.SetToCombination(1., fInternal_n, -1., farray);

  // checks on initial guess for deviatoric effective stress and damage
  if (fInternal[kDEQPe] <= 0.0)   fInternal[kDEQPe] = 1.e-10;
  if (fInternal[kALPH]  <= 0.0)   fInternal[kALPH] = 1.e-10;
  if (fInternal[kKAPP]  <= 0.0)   fInternal[kKAPP] = 1.e-10;
  if (fInternal[kDAMG]  <= 0.0)   fInternal[kDAMG] = DMG_LWB_LIMIT;
  if (fInternal[kDAMG]  >= 0.999) fInternal[kDAMG] = 0.999;
}

bool BCJHypoIsoDamageYC3D::IsSolnVariableNegative()
{
  // test for negative values of internal variables
  // note: volumetric strain DEQPh can be negative
  int i = 0;
  bool isNegative = false;
  while (i < fNumInternal && !isNegative)
    {
      if (i == kDEQPh) break;
      isNegative = (fInternal[i] < 0.0);
      i++;
    }

  // check for zero values of internal variables
  if (!isNegative) {
     if (fInternal[kDEQPe] <= 1.e-10)       fInternal[kDEQPe] = 1.e-10;
     if (fabs(fInternal[kDEQPh]) <= 1.e-10) fInternal[kDEQPh] = 1.e-10;
     if (fInternal[kALPH] <= 1.e-10)        fInternal[kALPH] = 1.e-10;
     if (fInternal[kKAPP] <= 1.e-10)        fInternal[kKAPP] = 1.e-10;
     if (fInternal[kDAMG] <= DMG_LWB_LIMIT) fInternal[kDAMG] = DMG_LWB_LIMIT;
     if (fInternal[kDAMG] >= 0.999)         fInternal[kDAMG] = 0.999;
  }

  return isNegative;
}

/* PRIVATE MEMBER FUNCTIONS */

void BCJHypoIsoDamageYC3D::ComputeInternalQntsRHS(const dArrayT& array)
{
  // equivalent stresses of porous material
  // ... coefficients A1 & A2 from void growth model
  fVoidGrowthModel->ACoefficients(array[kDAMG], fA1Dmg, fA2Dmg);
  fA1iDmg = 1./fA1Dmg;
  fA2iDmg = 1./fA2Dmg;

  // ... equivalent plastic strain increment of porous material
  fDEQP = sqrt(fA1iDmg*array[kDEQPe]*array[kDEQPe] + fA2iDmg*array[kDEQPh]*array[kDEQPh]);
  if (fabs(fDEQP) <= 1.e-6) fDEQP = 1.e-6;
 
  // ... equivalent stress of porous material
  fEQXi = fKineticEqn->h(fDEQP/fdt, array[kKAPP]);
  
  // ... deviatoric and hydrostatic equivalent stresses
  fEQValues[kEQXie] = fA1iDmg * array[kDEQPe]/fDEQP * fEQXi;
  fEQValues[kEQXih] = fA2iDmg * array[kDEQPh]/fDEQP * fEQXi;
 
  // no damage growth for possitive pressure (needs fixing)
  //if (fEQValues[kEQXih] < 0.e0) {
  //   fA2iDmg = 0.0;
  //   fDEQP = sqrt(fA1iDmg)*array[kDEQPe];
  //   if (fabs(fDEQP) <= 1.e-6) fDEQP = 1.e-6;
  //   fEQXi = fKineticEqn->h(fDEQP/fdt, array[kKAPP]);
  //   fEQValues[kEQXie] = fA1iDmg * array[kDEQPe]/fDEQP * fEQXi;
  //   fEQValues[kEQXih] = fEQXihTr;
  //} 

  // factor eta
  fEta = fMatProp[0]*array[kALPH]*array[kDEQPe] + fMatProp[2]*array[kALPH]*fdt;
  fEta = fEta/(1.+fEta);

  // "new" trial equivalente stress and unit normal in stress space
  fXiTr.SetToCombination(1., fSigTrDev, -2./3.*(1.-fEta), falpha_n);
  fEQXieTr = sqrt32*sqrt(fXiTr.ScalarProduct());
  fUnitNorm.SetToScaled(sqrt32/fEQXieTr, fXiTr);

  // factor (tensor) in integrated eqn for backstress: alpha=(1-eta)*X
  fX.SetToCombination(1., falpha_n, sqrt32*fMatProp[1]*array[kDEQPe], fUnitNorm);
  fXMag = sqrt(fX.ScalarProduct());
  if (fXMag <= 1.e-16) fXMag = 1.;  // to avoid NaN in fUnitM in first incr/iter
}

void BCJHypoIsoDamageYC3D::ComputeInternalQntsLHS(const dArrayT& array)
{
  // preliminary computations
  double betaStar = fKineticEqn->DhDeqpdot(fDEQP/fdt, array[kKAPP]) / fdt;
  double coefBeta = (betaStar*fDEQP-fEQXi)/fDEQP;
  double ratioE   = fEQValues[kEQXie]/fEQXi;
  double ratioH   = fEQValues[kEQXih]/fEQXi;

  // derivatives dA1/d() & dA2/d() from void growth model
  fVoidGrowthModel->ADerivCoefficients(array[kDAMG], fdA1Dmg, fdA2Dmg);

  // no damage growth for possitive pressure (needs fixing)
  //if (fA2iDmg == 0.0) { ratioH  = 0.0; fdA1Dmg = 0.0; fdA2Dmg = 0.0; }

  // derivative of (total) equivalent plastic strain increment wrt damage
  double dDEQPdDAMG = -0.5*fDEQP*(ratioE*ratioE*fdA1Dmg + ratioH*ratioH*fdA2Dmg);

  // derivatives of deviatoric equivalent stress wrt primary variables
  fdEQXie[0] = coefBeta*ratioE*ratioE + fA1iDmg*fEQXi/fDEQP;
  fdEQXie[1] = coefBeta*ratioE*ratioH;
  fdEQXie[2] = 0.;
  fdEQXie[3] = ratioE * fKineticEqn->DhDs(fDEQP/fdt, array[kKAPP]);
  fdEQXie[4] = coefBeta*ratioE*dDEQPdDAMG - fEQValues[kEQXie]*fA1iDmg*fdA1Dmg;

  // derivatives of hydrostatic equivalent stress wrt primary variables
  fdEQXih[0] = coefBeta*ratioE*ratioH;
  fdEQXih[1] = coefBeta*ratioH*ratioH + fA2iDmg*fEQXi/fDEQP;
  fdEQXih[2] = 0.;
  fdEQXih[3] = ratioH * fKineticEqn->DhDs(fDEQP/fdt, array[kKAPP]);
  fdEQXih[4] = coefBeta*ratioH*dDEQPdDAMG - fEQValues[kEQXih]*fA2iDmg*fdA2Dmg;

  // derivatives of factor Eta wrt primary variables
  fDEtaDDEQP = (1.-fEta)*(1.-fEta)*(fMatProp[0]*array[kALPH]);
  fDEtaDALPH = (1.-fEta)*(1.-fEta)*(fMatProp[0]*array[kDEQPe]+fMatProp[2]*fdt);
  // ...  other derivatives are zero

  // unit normal M: M=X/||X||
  fUnitM.SetToScaled(1./fXMag, fX);

  // derivative factor (tensor) for unit normal: dN/d()=A*dEta/d()
  falpha_n.ToMatrix(fmatx1);
  fUnitNorm.ToMatrix(fmatx2);
  fsymmatx1.SetToCombination(1., falpha_n, -dArrayT::Dot(fmatx1,fmatx2), fUnitNorm);
  fA.SetToScaled(1./(sqrt32*fEQXieTr), fsymmatx1);
}
