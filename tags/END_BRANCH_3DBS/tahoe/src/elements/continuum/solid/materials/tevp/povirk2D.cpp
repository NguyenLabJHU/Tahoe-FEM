/* $Id: povirk2D.cpp,v 1.9 2003-11-21 22:46:58 paklein Exp $ */
/* Created:  Harold Park (09/10/2001) */
#include "povirk2D.h"

#include <iostream.h>
#include <math.h>
#include "ifstreamT.h"
#include "FSMatSupportT.h"
#include "ElementCardT.h"

using namespace Tahoe;

const int kNumStateVariables = 4;
const int kNumOutput = 4;   // # of internal variables
const double kYieldTol = 1.0e-16;   // Yield stress criteria
const int kVoigt = 4;    // 4 components in 2D voigt notation
const int kNSD = 3;      // 3 2D stress components (11, 12=21, 22)
static const char* Labels[kNumOutput] = {
  "Temp",       // Temperature
  "Eff._Strain",    // effective strain
  "Eff._Stress",
  "Eff._Strain_Rate"};   // effective stress

/* constructor */
povirk2D::povirk2D(ifstreamT& in, const FSMatSupportT& support):
  FSSolidMatT(in, support),
  IsotropicT(in),
  Material2DT(in),        // Currently reads in plane strain from file...
  /* initialize references */
//  fRunState(ContinuumElement().RunState()),
//  fDt(ContinuumElement().ElementSupport().TimeStep()),
  fStress(2),
  fModulus(kVoigt),
  fLocVel(support.LocalArray(LocalArrayT::kVel)),
// fLocDisp(element.Displacements()),

  /* initialize work matrices */
  fGradV_2D(2),
  fGradV(3),
  fEbtot(0.0),
  fXxii(0.0),
  fEcc(0.0),
  fCtcon(0.0),
  fSpin(0.0),
  fTemperature(0.0),
  fSb(0.0),
  fEb(0.0),
  fEbRate(0.0),
  fStrainEnergyDensity(0.0),
  fCriticalStrain(0),
  fStill3D(3),
  fFtot_2D(2),
  fFtot(3),
  fDtot(3),          
  fF_temp(3),
  fPP(kVoigt),
  fDmat(kVoigt),
  fEP_tan(kVoigt),
  fStressMatrix(3),
  fStress3D(3),
  fSymStress2D(2),
  fStressArray(kVoigt),
  fSmlp(kVoigt),
  fVisc(0.0),
  fJ(0.0)
  
{

  /* initialize material constants */
  El_E = 2.0E11;
  El_V = .30;
  El_K = El_E / (3.0 * (1.0 - 2.0 * El_V));
  El_G = El_E / (2.0 * (1.0 + El_V));
  Sb0 = 2.0E9;
  Rho0 = 7830.0;
  Eb0 = Sb0 / El_E;
  Eb0tot = .001;
  BigN = .01;                // strain hardening exponent
  Smm = 70.0;                // rate sensitivity parameter

  /* initialize temperature parameters */
  Temp_0 = 293.0;
  Alpha_T = 11.2E-6;
  Delta = .8;
  Theta = .5;
  Kappa = 500.0;  
  Cp = 448.0;
  Chi = .9;
  Pcp = 3.0 * El_G;

  /* initialize damage parameters */
  Epsilon_1 = 4.0 * Eb0;
  Epsilon_2 = .3;
  Epsilon_rate = 4.0E4;
  Gamma_d = .002;
  SigCr = 6.0 * Sb0;

  /* used in temperature update */
  Xi = 1.0 / (Rho0 * Cp);
}

/* allocate element storage */
void povirk2D::PointInitialize(void)
{
  /* first ip only */
  if (CurrIP() == 0) AllocateElement(CurrentElement());
  LoadData(CurrentElement(), CurrIP());
  fInternal[kTemp] = 293.0;
}

/* required parameter flags */
bool povirk2D::NeedVel(void) const { return true; }

/* update internal variables */
void povirk2D::UpdateHistory(void)
{
  /* update if plastic */
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void povirk2D::ResetHistory(void)
{
  /* reset if plastic */
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void povirk2D::Print(ostream& out) const
{
  /* inherited */
  FSSolidMatT::Print(out);
  IsotropicT::Print(out);
  Material2DT::Print(out);
}

void povirk2D::PrintName(ostream& out) const
{
  /* inherited */
  FSSolidMatT::PrintName(out);
  out << "    Thermo-Elasto-Viscoplastic without fluid/damage model\n";
}

/* spatial description */
/* modulus */
const dMatrixT& povirk2D::C_IJKL(void)
{
  /* implement modulus here */
  return fModulus;   // Dummy - spatial description not desired
}

/* stress */
const dSymMatrixT& povirk2D::S_IJ(void)
{
  /* implement stress here */
  return fStress;    // Dummy - PK2 stress not desired
}

/* material description */
/* modulus */
const dMatrixT& povirk2D::c_ijkl(void)
{
  /* this function calculates the tangent modulus */

  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
  LoadData(element, ip);
  // Temporary work space arrays and matrices
  dMatrixT ata(kVoigt), lta(kVoigt);
  dArrayT diagU(kVoigt);

  ComputeDmat();   // Gives original elastic coefficient tensor
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  
  ComputePP();
  const double sb = fInternal[kSb];
  ComputeEbtotCtconXxii();

  /* calculate the tangent modulus - based on Pierce, 1987 */
  ata.Outer(fPP, fPP);
  lta.Outer(fPP, diagU);

  double k = fCtcon * 3.0 * El_K * Alpha_T * Chi * Xi * sb;
  fModulus = fDmat;
  fModulus.AddCombination(-fCtcon, ata, -k, lta);
  return fModulus;
}

/* stress */
const dSymMatrixT& povirk2D::s_ij(void)
{
  if (MaterialSupport().RunState() == GlobalT::kFormRHS)
  {
  	/* current time step */
	fDt = MaterialSupport().TimeStep();
  
    /* implement stress here - work with Kirchoff stress, then convert back
     * to Cauchy when necessary */
    /* Allocate all state variable space on first timestep, ie time = 0 */
    int ip = CurrIP();
    ElementCardT& element = CurrentElement();

    /* load data */
    LoadData(element, ip);
    iArrayT& flags = element.IntegerData();
    
    int checkplastic = CheckIfPlastic(element, ip);
    /* Incremental stress update part - if critical strain criteria not
     * exceeded */
    ComputeF();
    ComputeD();
    dArrayT sig_jrate(kVoigt), dtot(kVoigt), sts_dot(kVoigt);
    dMatrixT kirchoff_last(3);
    kirchoff_last = ArrayToMatrix(fTempKirchoff);
    const double spin = ComputeSpin();
    sig_jrate = 0.0;
    c_ijkl();               // Need the tangent modulus
    /* Flatten out the Rate of Deformation tensor into Voigt notation */    
    dtot[0] = fDtot(0,0);
    dtot[1] = fDtot(1,1);
    dtot[2] = 0.0;
    dtot[3] = fDtot(0,1);   // D is symmetric - could take (1,0)
    fModulus.Multx(dtot, sig_jrate);
    /* Check if plasticity has occurred yet */
    if (checkplastic == kIsPlastic) 
      sig_jrate -= ComputeEP_tan();
    
    /* Add the objective part */
    sts_dot[0] = sig_jrate[0] + 2.0 * kirchoff_last(0,1) * spin;
    sts_dot[1] = sig_jrate[1] - 2.0 * kirchoff_last(0,1) * spin;
    sts_dot[2] = sig_jrate[2];
    sts_dot[3] = sig_jrate[3] - spin * (kirchoff_last(0,0) - kirchoff_last(1,1));
    
    kirchoff_last(0,0) += sts_dot[0] * fDt;
    kirchoff_last(0,1) += sts_dot[3] * fDt;
    kirchoff_last(1,0) = kirchoff_last(0,1);
    kirchoff_last(1,1) += sts_dot[1] * fDt;
    kirchoff_last(2,2) += sts_dot[2] * fDt;
    fStress3D = Return3DStress(kirchoff_last);
    fTempKirchoff = MatrixToArray(fStress3D);    // This is the Kirchoff stress
    
    fStress3D /= fJ;
    fInternal[kSb] = ComputeEffectiveStress();
    fInternal[kTemp] = ComputeViscoTemperature();
    fInternal[kEb] = ComputeViscoEffectiveStrain();
    
    fStress.ReduceFrom3D(fStress3D);     // Take only 2D stress components
    // STORE CAUCHY STRESS HERE (2D version)
    fTempCauchy = fStress;
 
     /* store pressure */
	fInternal[kPressure] = fStress3D.Trace()/3.0;
 }
  else
  {
    LoadData(CurrentElement(), CurrIP());
    /* Extract cauchy stress from solution stage 6 */
    fStress = ArrayToSymMatrix2D(fTempCauchy);
  }
  /* return the stress */
  return fStress;
}

/* returns the strain energy density for the specified strain */
double povirk2D::StrainEnergyDensity(void)
{
  /* compute strain energy density here */
  return fStrainEnergyDensity;
}

int povirk2D::NumOutputVariables(void) const { return kNumOutput; }
void povirk2D::OutputLabels(ArrayT<StringT>& labels) const
{
  /* set size */
  labels.Dimension(kNumOutput);

  /* copy labels - WHY? */
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void povirk2D::ComputeOutput(dArrayT& output)
{
  /* Currently assuming that UpdateHistory is called before ComputeOutput */
  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
  LoadData(element, ip);
  output[0] = fInternal[kTemp];        // Temperature
  output[1] = fInternal[kEb];          // Effective strain
  output[2] = fInternal[kSb];          // Effective stress
  output[3] = fEbRate;                 // Effective strain rate
}

/*******************************************************************
Computational routines start here
*******************************************************************/

void povirk2D::ComputeF(void)
{
  /* Compute deformation gradient, Jacobian, inverse of deformation gradient */
  fFtot = fF_temp = 0.0;
  dMatrixT* tempf = &fFtot;
  dMatrixT* f_inv = &fF_temp;
  double* j = &fJ;
  fFtot_2D = F();
  (*tempf).Rank2ExpandFrom2D(fFtot_2D);
  (*tempf)(2,2) = 1.0;

  *j = fFtot.Det();
  (*f_inv).Inverse(fFtot);
}

void povirk2D::ComputeD(void)
{
	if (!fLocVel) throw ExceptionT::kGeneralFail;

	/* Compute rate of deformation */
	fDtot = 0.0;
	dMatrixT* tempd = &fDtot;
	if (!FSMatSupport().ComputeGradient_reference(*fLocVel, fGradV_2D)) throw ExceptionT::kGeneralFail;
	fGradV.Rank2ExpandFrom2D(fGradV_2D);
	(*tempd).MultAB(fGradV, fF_temp, 0);
	(*tempd).Symmetrize();
}

double povirk2D::ComputeSpin(void)
{
	if (!fLocVel) throw ExceptionT::kGeneralFail;

	/* Compute the spin scalar */
	fSpin = 0.0;
	if (!FSMatSupport().ComputeGradient_reference(*fLocVel, fGradV_2D)) throw ExceptionT::kGeneralFail;
	fGradV.Rank2ExpandFrom2D(fGradV_2D);
	fSpin = fGradV(0,0) * fF_temp(0,1) + fGradV(0,1) * fF_temp(1,1);
	fSpin = fSpin - fGradV(1,0) * fF_temp(0,0) - fGradV(1,1) * fF_temp(1,0);
	fSpin *= .5;
	return fSpin;
}

double povirk2D::ComputeViscoTemperature(void)
{
  /* Case where fluid model was not used - viscoplasticity */
  const double sb = fInternal[kSb];
  const double temp_last = fInternal[kTemp];
  double temp_rate = Chi * Xi * fEbtot * sb;
  fTemperature = temp_rate * fDt + temp_last;

  return fTemperature;
}

double povirk2D::ComputeViscoEffectiveStrain(void)
{
  /* If fluid model / critical strain is not used */
  const double eb_last = fInternal[kEb];
  double ecc = ComputeEcc();
  
  /* access necessary data */
  double *ebtot_c = &fEbRate;
  *ebtot_c = fEbtot / (1.0 + fXxii) + fCtcon * ecc;
  fEb = eb_last + fDt * (*ebtot_c);

  return fEb;
}

double povirk2D::ComputeEffectiveStress(void)
{
  /* Computes the effective stress - whether damage/fluid model is used
   * is apparently not relevant - it's computed the same way */
  dMatrixT temp_stress(3);
  temp_stress = ArrayToMatrix(fTempKirchoff);
  double trace = temp_stress.Trace() / 3.0;
  double temp1 = pow(temp_stress(0,0) - trace, 2);
  double temp2 = pow(temp_stress(1,1) - trace, 2);
  double temp3 = pow(temp_stress(2,2) - trace, 2);
  fSb = 1.5 * (temp1 + temp2 + temp3) + 3.0 * pow(temp_stress(0,1), 2);
  fSb = sqrt(fSb);

  return fSb;
}

void povirk2D::ComputeEbtotCtconXxii(void)
{
  /* compute Ctcon - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  double* ctcon = &fCtcon;
  double* xxii = &fXxii;
  double* ebtot = &fEbtot;
  const double sb = fInternal[kSb];
                              
  if (sb <= kYieldTol)
  {
    *ctcon = 0.0;
    *ebtot = 0.0;
    *xxii = 0.0;
  }
  else
  {
    const double eb = fInternal[kEb];	
    const double temp = fInternal[kTemp];	
    double pCp2 = dArrayT::Dot(fSmlp, fPP);
    pCp2 += Pcp;
    double gsoft = Sb0 * pow(1.0 + eb/Eb0, BigN);
    gsoft *= (1.0 - Delta * (exp((temp - Temp_0)/Kappa) - 1.0));
    double reg = sb / gsoft;
    *ebtot = Eb0tot * pow(reg, Smm);
    double pE_ptau = Smm * (*ebtot) / sb;
    double dG_Eb = BigN * gsoft / (eb + Eb0);
    double dG_T = -(Sb0 * Delta / Kappa) * (pow(1.0 + eb / Eb0, BigN) * exp((temp - Temp_0) / Kappa));
    double pE_evt = -dG_Eb * (sb / gsoft);
    double pE_Tvt = -dG_T * (sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * Chi * Xi * sb;

    *xxii = Theta * fDt * hh * pE_ptau;
    *ctcon = *xxii / ((1.0 + (*xxii)) * hh);
  }
}

void povirk2D::ComputeSmlp(void)
{
  /* used in ComputePP */
  /* compute the deviatoric Kirchoff stress */
  fSmlp = 0.0;
  dArrayT* smlp = &fSmlp;
  const double sb = fInternal[kSb];
  if (sb <= kYieldTol) 
    *smlp = 0.0;
  else
  {
    dMatrixT temp_stress(3); 
    temp_stress = ArrayToMatrix(fTempKirchoff);
    double trace_KH = temp_stress.Trace() / 3.0;
    (*smlp)[0] = temp_stress(0,0) - trace_KH;
    (*smlp)[1] = temp_stress(1,1) - trace_KH;
    (*smlp)[2] = temp_stress(2,2) - trace_KH;
    (*smlp)[3] = temp_stress(0,1);   // Stored in Voigt notation
    *smlp *= 1.5 / sb;
  }
}

void povirk2D::ComputePP(void)
{
  fPP = 0.0;
  dArrayT* pp = &fPP;
  ComputeSmlp();
  fDmat.Multx(fSmlp, *pp);
}

double povirk2D::ComputeEcc(void)
{
  /* Access ecc */
  fEcc = 0.0;
  if (fInternal[kSb] <= kYieldTol)   // If hasn't yielded yet...
    fEcc = 0.0;
  else
  {
    for (int i = 0; i < 3; i++) 
      fEcc += fPP[i] * fDtot(i,i);

    fEcc += fPP[3] * fDtot(0,1);
  }
  return fEcc;
}

void povirk2D::ComputeDmat(void)
{
  /* computes the original elastic coefficient tensor */
  fDmat = 0.0;
  dMatrixT* dmat = &fDmat;
  double ed = El_E * (1.0 - El_V) / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double es = El_E * El_V / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double g0 = El_E / (1.0 + El_V);
  (*dmat)(0,0) = (*dmat)(1,1) = (*dmat)(2,2) = ed;
  (*dmat)(0,1) = (*dmat)(1,0) = (*dmat)(2,0) = (*dmat)(2,1) = es;
  (*dmat)(0,2) = (*dmat)(0,3) = (*dmat)(1,2) = (*dmat)(1,3) = 0.0;
  (*dmat)(2,3) = (*dmat)(3,0) = (*dmat)(3,1) = (*dmat)(3,2) = 0.0;
  (*dmat)(3,3) = g0;

}

dArrayT& povirk2D::ComputeEP_tan(void)
{
  /* computes the modulus correction if plasticity has occurred */
  fEP_tan = 0.0;
  dArrayT diagU(kVoigt);
  const double sb = fInternal[kSb];
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  

  for (int i = 0; i < 4; i++) {
    /* EP_tan is the plastic corrector to the tangent modulus */
    fEP_tan[i] = (fEbtot / (1.0 + fXxii)) * (fPP[i] + 3.0 * El_K * Alpha_T * Xi * Chi * sb * diagU[i]);  
  }

  return fEP_tan;
}

int povirk2D::CheckIfPlastic(ElementCardT& element, int ip)
{
  /* Checks to see if the gauss point has gone plastic yet via a
   * test on the effective stress */
  /* plastic */
  iArrayT& flags = element.IntegerData();
  if (fInternal[kSb] > kYieldTol)
  {
    flags[ip] = kIsPlastic;    // Has gone plastic
    return 0;
  }
  /* elastic */
  else
  {
    flags[ip] = kIsElastic;   // Hasn't gone plastic yet
    return 1;
  }
}

void povirk2D::AllocateElement(ElementCardT& element)
{
  /* return a pointer to a new plastic element object constructed with
   * the data from element */
  /* determine storage */
  int i_size = 0;
  int d_size = 0;
  i_size += 2 * fNumIP;              // 2 flags per IP:  critical strain
                                      // and check for plasticity
  d_size += kNumStateVariables * fNumIP;     // 3 internal variables to track
  d_size += kVoigt * fNumIP;         // 4 non-zero stress components:
                                      // Sig11, Sig12=Sig21, Sig22 and Sig33
  d_size += kNSD * fNumIP;           // 3 2D symmetric components (Sig11, Sig12, Sig22)
  /* construct new plastic element */
  element.Dimension(i_size, d_size);

  /* first set of flags for plasticity criterion */
  for (int ip = 0; ip < fNumIP; ip++)
    (element.IntegerData())[ip] = kIsElastic;
  
  /* second set of flags for critical strain / model switch criterion */
  for (int ip = fNumIP; ip < 2 * fNumIP; ip++)
    (element.IntegerData())[ip] = kTevp;

  element.DoubleData() = 0.0;
}

void povirk2D::LoadData(const ElementCardT& element, int ip)
{
  /* load element data for the specified integration point */
  /* check */
  if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

  int dex = ip * kVoigt;     // 4 non-zero stress components (11, 12, 22, 33)
  int dex2 = ip * kNSD;      // 3 non-zero 2D stress components (11, 12=21, 22)
  int offset = fNumIP * 4;
  int offset2 = kNSD * offset / kVoigt;
  /* fetch arrays */
  const dArrayT& d_array = element.DoubleData();
  fTempKirchoff.Alias(kVoigt, &d_array[dex]);
  fTempCauchy.Alias(kNSD, &d_array[offset + dex2]);
  fInternal.Alias(kNumStateVariables, &d_array[offset + offset2 + ip * kNumStateVariables]); 
}

void povirk2D::Update(ElementCardT& element)
{
  /* get flags */
  iArrayT& flags = element.IntegerData();
  /* check if reset state (is same for all ip) */
  if (flags[0] == kReset)
  {
      flags = kIsElastic;          // don't update again
      return;
  }
}

void povirk2D::Reset(ElementCardT& element)
{
  /* resets to the last converged solution */
  /* flag to not update again */
  (element.IntegerData()) = kReset;
}

dArrayT& povirk2D::MatrixToArray(const dSymMatrixT& StressMatrix)
{
  /* Flattens Kirchoff stress matrix into array form for internal variable
   * storage */
  fStressArray[0] = StressMatrix[0];        // Sigma 11
  fStressArray[3] = StressMatrix[5];        // Sigma 12
  fStressArray[1] = StressMatrix[1];        // Sigma 22
  fStressArray[2] = StressMatrix[2];        // Sigma 33
  
  return fStressArray;
}

dMatrixT& povirk2D::ArrayToMatrix(const dArrayT& StressArray)
{
  /* Expands internal variable stress array to matrix form */
  fStressMatrix = 0.0;
  fStressMatrix(0,0) = StressArray[0];       // Sigma 11
  fStressMatrix(1,0) = fStressMatrix(0,1) = StressArray[3];   // Sigma 12
  fStressMatrix(1,1) = StressArray[1];       // Sigma 22
  fStressMatrix(2,2) = StressArray[2];       // Sigma 33
  return fStressMatrix;
}

dSymMatrixT& povirk2D::Return3DStress(const dMatrixT& StressMatrix)
{
  /* Takes 3D matrix and converts to 3D symmetric matrix - necessary
   * because canned functions depend on fNumSD to convert */
  fStill3D[0] = StressMatrix(0,0);      // Sigma 11
  fStill3D[1] = StressMatrix(1,1);      // Sigma 22
  fStill3D[2] = StressMatrix(2,2);      // Sigma 33
  fStill3D[3] = StressMatrix(1,2);      // Sigma 23
  fStill3D[4] = StressMatrix(0,2);      // Sigma 13
  fStill3D[5] = StressMatrix(0,1);      // Sigma 12
  return fStill3D;
}

dSymMatrixT& povirk2D::ArrayToSymMatrix2D(const dArrayT& StressArray)
{
  /* Because dSymMatrixT is derived from dArrayT, cannot assign fStress = StressArray */
  fSymStress2D[0] = StressArray[0];     // Sigma 11
  fSymStress2D[1] = StressArray[1];     // Sigma 22
  fSymStress2D[2] = StressArray[2];     // Sigma 33

  return fSymStress2D;
}

