/* $Id: tevp3D.cpp,v 1.16 2003-11-21 22:46:58 paklein Exp $ */
/* created:  Harold Park (06/25/2001) */
#include "tevp3D.h"
#include <iostream.h>
#include <math.h>
#include "ifstreamT.h"
#include "FSMatSupportT.h"
#include "ElementCardT.h"

using namespace Tahoe;

/* element output data */
const int kNumOutput = 3;   // # of internal variables
const double kYieldTol = 1.0e-16;   // Yield stress criteria
const int kVoigt = 6;    // 6 stress components in 3D voigt notation
static const char* Labels[kNumOutput] = {
  "Temp",       // Temperature
  "Eff._Strain",    // effective strain
  "Eff._Stress"};   // effective stress

/* constructor */
tevp3D::tevp3D(ifstreamT& in, const FSMatSupportT& support):
  FSSolidMatT(in, support),
  IsotropicT(in),
  /* initialize references */
//  fRunState(ContinuumElement().RunState()),
//  fDt(ContinuumElement().ElementSupport().TimeStep()),
  fStress(3),
  fModulus(kVoigt),
  fLocVel(support.LocalArray(LocalArrayT::kVel)),
//  fLocDisp(element.Displacements()),

  /* initialize work matrices */
  fGradV(3),
  fEbtot(0.0),
  fXxii(0.0),
  fEcc(0.0),
  fCtcon(0.0),
  fSpin(3),
  fTemperature(0.0),
  fSb(0.0),
  fEb(0.0),
  fStrainEnergyDensity(0.0),
  fCriticalStrain(0),
  fStill3D(3),
  fFtot(3),
  fDtot(3),          
  fF_temp(3),
  fPP(kVoigt),
  fDmat(kVoigt),
  fEP_tan(kVoigt),
  fStressMatrix(3),
  fStress3D(3),
  fSymStress3D(3),
  fStressArray(kVoigt),
  fSmlp(kVoigt),
  fVisc(0.0),
  fJ(0.0)
  
{
	/* read parameters from input stream */
	in >> Mu_d;
	fVisc = Mu_d;

	/* check values */
	if (Mu_d < 0.0) 
	{
		cout << "\n tevp3D::tevp3D: Mu must be > 0: " << Mu_d << endl;
		throw ExceptionT::kBadInputValue;
	}

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
  //Mu_d = 500.0;
  SigCr = 6.0 * Sb0;

  /* used in temperature update */
  Xi = 1.0 / (Rho0 * Cp);
}

/* allocate element storage */
void tevp3D::PointInitialize(void)
{
  /* first ip only */
  if (CurrIP() == 0) AllocateElement(CurrentElement());
  LoadData(CurrentElement(), CurrIP());
  fInternal[kTemp] = 293.0;
}

/* required parameter flags */
bool tevp3D::NeedVel(void) const { return true; }

/* update internal variables */
void tevp3D::UpdateHistory(void)
{
  /* update if plastic */
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void tevp3D::ResetHistory(void)
{
  /* reset if plastic */
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void tevp3D::Print(ostream& out) const
{
  /* inherited */
  FSSolidMatT::Print(out);
  IsotropicT::Print(out);
}

void tevp3D::PrintName(ostream& out) const
{
  /* inherited */
  FSSolidMatT::PrintName(out);
  out << "    3D Thermo-Elasto-Viscoplastic\n";
}

/* spatial description */
/* modulus */
const dMatrixT& tevp3D::C_IJKL(void)
{
  /* implement modulus here */
  return fModulus;   // Dummy - spatial description not desired
}

/* stress */
const dSymMatrixT& tevp3D::S_IJ(void)
{
	cout << "\n tevp3D::S_IJ: not implemented" << endl;
	throw ExceptionT::kGeneralFail;

	/* implement stress here */
	return fStress;
}

/* material description */
/* modulus */
const dMatrixT& tevp3D::c_ijkl(void)
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
  diagU[3] = diagU[4] = diagU[5] = 0.0;  
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
const dSymMatrixT& tevp3D::s_ij(void)
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
 
    if (flags[ip + fNumIP] == kFluid) {
      /* Fluid model part - if critical strain criteria is exceeded */
      ComputeF();
      ComputeD();
      const double temp = fInternal[kTemp];   // Use the PREVIOUS temperature
      double cm = -Gamma_d * El_E * (1.0 - fJ + Alpha_T * (temp - Temp_0));
      cm /= (fJ * (1.0 - El_V));
      
      /* add dynamic viscosity */
      double *mu_d = &Mu_d;
      *mu_d = fVisc * exp((temp - Temp_0) / Temp_0);
      dMatrixT eye_cm(3), dtemp(3);
      eye_cm = 0.0;
      eye_cm.PlusIdentity(1.0);
      eye_cm *= cm;
      fDtot.ToMatrix(dtemp);          
      dtemp *= Mu_d;
      eye_cm += dtemp;
      fStress = Return3DStress(eye_cm);
      fTempKirchoff = MatrixToArray(fStress);
      fStress /= fJ;           // Return the Cauchy, NOT Kirchoff stress!!!
      fInternal[kSb] = ComputeEffectiveStress();
      fInternal[kTemp] = ComputeFluidTemperature();
      fInternal[kEb] = ComputeFluidEffectiveStrain();
    }  
    else {
      /* Incremental stress update part - if critical strain criteria not
       * exceeded */
      ComputeF();
      ComputeD();
      dArrayT sig_jrate(kVoigt);
      dMatrixT kirchoff_last(3), spin(3), sts_dot(3), part1(3), part2(3);
      kirchoff_last = ArrayToMatrix(fTempKirchoff);
      spin = ComputeSpin();
      sig_jrate = 0.0;
      c_ijkl();               // Need the tangent modulus

      fModulus.Multx(fDtot, sig_jrate);
      /* Check if plasticity has occurred yet */
      if (checkplastic == kIsPlastic) 
        sig_jrate -= ComputeEP_tan();
      
      /* Add the objective part */
      sts_dot = ArrayToMatrix(sig_jrate);
      part1.MultAB(spin, kirchoff_last);
      part2.MultAB(kirchoff_last, spin);
      sts_dot += part1;
      sts_dot -= part2;
      sts_dot *= fDt;
      kirchoff_last += sts_dot;

      fStress = Return3DStress(kirchoff_last);
      fTempKirchoff = MatrixToArray(fStress);    // This is the Kirchoff stress
      
      fStress /= fJ;
      fInternal[kSb] = ComputeEffectiveStress();
      fInternal[kTemp] = ComputeViscoTemperature();
      fInternal[kEb] = ComputeViscoEffectiveStrain();
    }
    
    /* store cauchy stress here */
    fTempCauchy = fStress;
    CheckCriticalCriteria(element, ip);
  }
  else
  {
    LoadData(CurrentElement(), CurrIP());
    /* Extract cauchy stress from solution stage 6 */
    fStress = ArrayToSymMatrix3D(fTempCauchy);
  }
  /* return the stress */
  return fStress;
}

/* returns the strain energy density for the specified strain */
double tevp3D::StrainEnergyDensity(void)
{
  /* compute strain energy density here */
  return fStrainEnergyDensity;
}

int tevp3D::NumOutputVariables(void) const { return kNumOutput; }
void tevp3D::OutputLabels(ArrayT<StringT>& labels) const
{
  /* set size */
  labels.Dimension(kNumOutput);

  /* copy labels - WHY? */
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void tevp3D::ComputeOutput(dArrayT& output)
{
  /* Currently assuming that UpdateHistory is called before ComputeOutput */
  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
  LoadData(element, ip);
  output[0] = fInternal[kTemp];        // Temperature
  output[1] = fInternal[kEb];          // Effective strain
  output[2] = fInternal[kSb];          // Effective stress
}

/*******************************************************************
Computational routines start here
*******************************************************************/

void tevp3D::ComputeF(void)
{
  /* Compute deformation gradient, Jacobian, inverse of deformation gradient */
  fFtot = fF_temp = 0.0;
  dMatrixT* tempf = &fFtot;
  dMatrixT* f_inv = &fF_temp;
  double* j = &fJ;
  *tempf = F();

  *j = fFtot.Det();
  (*f_inv).Inverse(fFtot);
}

void tevp3D::ComputeD(void)
{
	if (!fLocVel) throw ExceptionT::kGeneralFail;

  /* Compute rate of deformation, put in in 3D symmetric stress array form */
  fDtot = 0.0;
  dSymMatrixT* smalld = &fDtot;
  dSymMatrixT tempd(3);
  dMatrixT yada(3);
  if (!FSMatSupport().ComputeGradient_reference(*fLocVel, fGradV)) throw ExceptionT::kGeneralFail;
  yada.MultAB(fGradV, fF_temp, 0);
  (*smalld) = tempd.Symmetrize(yada);
}

dMatrixT& tevp3D::ComputeSpin(void)
{
	if (!fLocVel) throw ExceptionT::kGeneralFail;

  /* Compute the spin tensor */
  fSpin = 0.0;
  if (!FSMatSupport().ComputeGradient_reference(*fLocVel, fGradV)) throw ExceptionT::kGeneralFail;
  dMatrixT yada(3);
  yada.MultAB(fGradV, fF_temp, 0);
  double temp1 = .5 * (yada(0,1) - yada(1,0));
  double temp2 = .5 * (yada(0,2) - yada(2,0));
  double temp3 = .5 * (yada(1,2) - yada(2,1));
  fSpin(0,1) = temp1;
  fSpin(0,2) = temp2;
  fSpin(1,2) = temp3;
  fSpin(1,0) = -temp1;
  fSpin(2,0) = -temp2;
  fSpin(2,1) = -temp3;

  return fSpin;
}

double tevp3D::ComputeFluidTemperature(void)
{
  /* Compute the output temperature - 2 different methods of computing, 
   * which depends upon whether fluid model was used or not */

  const double temp_last = fInternal[kTemp];
  
  /* Case where fluid model was used */
  dMatrixT temp_stress(3); 
  temp_stress = ArrayToMatrix(fTempKirchoff);
  double wpdot = temp_stress(0,0) * fDtot[0] + temp_stress(1,1) * fDtot[1] + temp_stress(2,2) * fDtot[2];
  wpdot += (2.0 * Mu_d * (fDtot[3] * fDtot[3] + fDtot[4] * fDtot[4] + fDtot[5] * fDtot[5]));
  double temp_rate = Chi * Xi * wpdot;
  fTemperature = temp_rate * fDt + temp_last;

  return fTemperature;
}

double tevp3D::ComputeViscoTemperature(void)
{
  /* Case where fluid model was not used - viscoplasticity */
  const double sb = fInternal[kSb];
  const double temp_last = fInternal[kTemp];
  double temp_rate = Chi * Xi * fEbtot * sb;
  fTemperature = temp_rate * fDt + temp_last;

  return fTemperature;
}

double tevp3D::ComputeFluidEffectiveStrain(void)
{
  /* Computes the effective strain - 2 different methods of computing,
   * which depends upon whether fluid model was used or not */
  const double eb_last = fInternal[kEb];
  double temp1 = fDtot[0] * fDtot[0] + fDtot[1] * fDtot[1] + fDtot[2] * fDtot[2];
  double temp2 = fDtot[3] * fDtot[3] + fDtot[4] * fDtot[4] + fDtot[5] * fDtot[5];
  double ebar = 1.5 * temp1 + 2.0 * temp2;
  ebar = sqrt(ebar);
  fEb = ebar * fDt + eb_last;

  return fEb;
}

double tevp3D::ComputeViscoEffectiveStrain(void)
{
  /* If fluid model / critical strain is not used */
  const double eb_last = fInternal[kEb];
  double ecc = ComputeEcc();
  
  /* access necessary data */
  double ebtot_c = fEbtot / (1.0 + fXxii) + fCtcon * ecc;
  fEb = eb_last + fDt * ebtot_c;

  return fEb;
}

double tevp3D::ComputeEffectiveStress(void)
{
  /* Computes the effective stress - whether damage/fluid model is used
   * is apparently not relevant - it's computed the same way */
  dMatrixT temp_stress(3);
  temp_stress = ArrayToMatrix(fTempKirchoff);
  double trace = temp_stress.Trace() / 3.0;
  double temp1 = pow(temp_stress(0,0) - trace, 2);
  double temp2 = pow(temp_stress(1,1) - trace, 2);
  double temp3 = pow(temp_stress(2,2) - trace, 2);
  fSb = 1.5 * (temp1 + temp2 + temp3) + 3.0 * (pow(temp_stress(0,1), 2) + pow(temp_stress(2,0), 2) + pow(temp_stress(1,2), 2));
  fSb = sqrt(fSb);

  return fSb;
}

void tevp3D::CheckCriticalCriteria(ElementCardT& element, int ip)
{
  /* Returns an indicator to determine whether critical strain criteria
   * has been met, and switch to fluid model happens next time step */
  iArrayT& flags = element.IntegerData();
  int* criticalstrain = &fCriticalStrain;

  /* if already fluid, no need to check criterion */
  if (flags[ip + fNumIP] == kFluid)
    *criticalstrain = 1;
  else if (flags[ip + fNumIP] == kCrack) 
    *criticalstrain = 2;
  else
  {  
    const double eb = fInternal[kEb];
    double ebar_cr = Epsilon_1 + (Epsilon_2 - Epsilon_1) * Epsilon_rate;
    ebar_cr /= (Epsilon_rate + fEbtot);
    if (eb >= ebar_cr)
    {
      flags[ip + fNumIP] = kFluid;
      *criticalstrain = 1;        // Indicator to switch to fluid model
    }
    else
    {
      flags[ip + fNumIP] = kTevp; 
      *criticalstrain = 0;
    }
  } 
}

void tevp3D::ComputeEbtotCtconXxii(void)
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

void tevp3D::ComputeSmlp(void)
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
    (*smlp)[3] = temp_stress(2,1);   // S23 - stored in Voigt notation 
    (*smlp)[4] = temp_stress(2,0);   // S31
    (*smlp)[5] = temp_stress(0,1);   // S12
    *smlp *= 1.5 / sb;
  }
}

void tevp3D::ComputePP(void)
{
  fPP = 0.0;
  dArrayT* pp = &fPP;
  ComputeSmlp();
  fDmat.Multx(fSmlp, *pp);
}

double tevp3D::ComputeEcc(void)
{
  /* Access ecc */
  fEcc = 0.0;
  if (fInternal[kSb] <= kYieldTol)   // If hasn't yielded yet...
    fEcc = 0.0;
  else
    fEcc = dArrayT::Dot(fPP, fDtot);

  return fEcc;
}

void tevp3D::ComputeDmat(void)
{
  /* computes the original elastic coefficient tensor */
  fDmat = 0.0;
  dMatrixT* dmat = &fDmat;
  double el_lambda = El_K - (2 / 3) * El_G;
  double C11 = el_lambda + 2 * El_G;
  double C44 = 2 * El_G;
  (*dmat)(0,0) = (*dmat)(1,1) = (*dmat)(2,2) = C11;
  (*dmat)(0,1) = (*dmat)(0,2) = (*dmat)(1,0) = (*dmat)(1,2) = (*dmat)(2,0) = (*dmat)(2,1) = el_lambda;
  (*dmat)(3,3) = (*dmat)(4,4) = (*dmat)(5,5) = C44;
}

dArrayT& tevp3D::ComputeEP_tan(void)
{
  /* computes the modulus correction if plasticity has occurred */
  fEP_tan = 0.0;
  dArrayT diagU(kVoigt);
  const double sb = fInternal[kSb];
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = diagU[4] = diagU[5] = 0.0;  

  for (int i = 0; i < kVoigt; i++) {
    /* EP_tan is the plastic corrector to the tangent modulus - based on Peirce, 1987*/
    fEP_tan[i] = (fEbtot / (1.0 + fXxii)) * (fPP[i] + 3.0 * El_K * Alpha_T * Xi * Chi * sb * diagU[i]);  
  }

  return fEP_tan;
}

int tevp3D::CheckIfPlastic(ElementCardT& element, int ip)
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

void tevp3D::AllocateElement(ElementCardT& element)
{
  /* return a pointer to a new plastic element object constructed with
   * the data from element */
  /* determine storage */
  int i_size = 0;
  int d_size = 0;
  i_size += 2 * fNumIP;              // 2 flags per IP:  critical strain
                                      // and check for plasticity
  d_size += kNumOutput * fNumIP;     // 3 internal variables to track
  d_size += kVoigt * fNumIP;         // 6 symmetric stress components:
                                      // S11, S22, S33, S23, S13, S21
  d_size += kVoigt * fNumIP;           // 6 3D symmetric components (S11, S22, S33, S23, S13, S12)
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

void tevp3D::LoadData(const ElementCardT& element, int ip)
{
  /* load element data for the specified integration point */
  /* check */
  if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

  int dex = ip * kVoigt;     // 6 non-zero 3Dstress components (11, 22, 33, 23, 13, 21)
  int offset = fNumIP * kVoigt;
  /* fetch arrays */
  const dArrayT& d_array = element.DoubleData();
  fTempKirchoff.Alias(kVoigt, &d_array[dex]);
  fTempCauchy.Alias(kVoigt, &d_array[offset + dex]);
  fInternal.Alias(kNumOutput, &d_array[2 * offset + ip * kNumOutput]);
}

void tevp3D::Update(ElementCardT& element)
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

void tevp3D::Reset(ElementCardT& element)
{
  /* resets to the last converged solution */
  /* flag to not update again */
  (element.IntegerData()) = kReset;
}

dArrayT& tevp3D::MatrixToArray(const dSymMatrixT& StressMatrix)
{
  /* Flattens Kirchoff stress matrix into array form for internal variable
   * storage */
  fStressArray[0] = StressMatrix[0];        // Sigma 11
  fStressArray[1] = StressMatrix[1];        // Sigma 22
  fStressArray[2] = StressMatrix[2];        // Sigma 33
  fStressArray[3] = StressMatrix[3];        // Sigma 23
  fStressArray[4] = StressMatrix[4];        // Sigma 13
  fStressArray[5] = StressMatrix[5];        // Sigma 12

  return fStressArray;
}

dMatrixT& tevp3D::ArrayToMatrix(const dArrayT& StressArray)
{
  /* Expands internal variable stress array to matrix form */
  fStressMatrix = 0.0;
  fStressMatrix(0,0) = StressArray[0];       // Sigma 11
  fStressMatrix(1,1) = StressArray[1];       // Sigma 22
  fStressMatrix(2,2) = StressArray[2];       // Sigma 33
  fStressMatrix(2,1) = fStressMatrix(1,2) = StressArray[3];   // Sigma 23
  fStressMatrix(0,2) = fStressMatrix(2,0) = StressArray[4];   // Sigma 13
  fStressMatrix(1,0) = fStressMatrix(0,1) = StressArray[5];   // Sigma 12

  return fStressMatrix;
}

dSymMatrixT& tevp3D::Return3DStress(const dMatrixT& StressMatrix)
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

dSymMatrixT& tevp3D::ArrayToSymMatrix3D(const dArrayT& StressArray)
{
  /* Because dSymMatrixT is derived from dArrayT, cannot assign fStress = StressArray */
  fSymStress3D[0] = StressArray[0];     // Sigma 11
  fSymStress3D[1] = StressArray[1];     // Sigma 22
  fSymStress3D[2] = StressArray[2];     // Sigma 33
  fSymStress3D[3] = StressArray[3];     // Sigma 23
  fSymStress3D[4] = StressArray[4];     // Sigma 13
  fSymStress3D[5] = StressArray[5];     // Sigma 12

  return fSymStress3D;
}

