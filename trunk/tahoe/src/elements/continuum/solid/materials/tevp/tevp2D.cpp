/* $Id: tevp2D.cpp,v 1.5 2001-05-03 02:48:20 hspark Exp $ */
/* Implementation file for thermo-elasto-viscoplastic material subroutine */
/* Created:  Harold Park (04/04/2001) */
/* Last Updated:  Harold Park (05/02/2001) */

#include "tevp2D.h"
#include <iostream.h>
#include <math.h>
#include "ElasticT.h"
#include "ShapeFunctionT.h"
#include "FEManagerT.h"
#include "ElementCardT.h"

/* element output data */
const int kNumOutput = 3;   // # of internal variables
const double kYieldTol = 1.0e-16;   // Yield stress criteria
const int kVoigt = 4;    // 4 components in 2D voigt notation
const int kNSD = 3;      // 3 2D stress components (11, 12=21, 22)
static const char* Labels[kNumOutput] = {
  "Temp",       // Temperature
  "Eff._Strain",    // effective strain
  "Eff._Stress"};   // effective stress

/* constructor */
tevp2D::tevp2D(ifstreamT& in, const ElasticT& element):
  FDStructMatT(in, element),
  IsotropicT(in),
  Material2DT(in),        // Currently reads in plane strain from file...
  /* initialize references */
  fRunState(ContinuumElement().RunState()),
  fDt(ContinuumElement().FEManager().TimeStep()),
  fTime(ContinuumElement().FEManager().Time()),
  fStress(2),
  fModulus(kVoigt),
  fShapes(element.ShapeFunction()),
  fLocVel(element.Velocities()),
  fLocDisp(element.Displacements()),

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
  fStrainEnergyDensity(0.0),
  fCriticalStrain(0),
  fStill3D(3),
  fKirchoff(3),
  fFtot_2D(2),
  fFtot(3),
  fDtot(3),          
  fF_temp(3),
  fFinv_2D(3),
  fPP(kVoigt),
  fDmat(kVoigt),
  fEP_tan(kVoigt),
  fStressMatrix(3),
  fStress3D(3),
  fSymStress2D(2),
  fStressArray(kVoigt),
  fSmlp(kVoigt)
  
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
  Mu_d = 500.0;
  SigCr = 6.0 * Sb0;

  /* used in temperature update */
  Xi = 1.0 / (Rho0 * Cp);
  cout << "CONSTRUCTOR INSTANTIATED AGAIN!!" << endl;
}

/* allocate element storage */
void tevp2D::PointInitialize(void)
{
  cout << "PointInitialize called" << endl;
  /* first ip only */
  if (CurrIP() == 0) AllocateElement(CurrentElement());
  LoadData(CurrentElement(), CurrIP());
  fInternal[kTemp] = 293.0;
}

/* required parameter flags */
bool tevp2D::NeedVel(void) const { return true; }

/* update internal variables */
void tevp2D::UpdateHistory(void)
{
  /* update if plastic */
  cout << "UpdateHistory called" << endl;
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void tevp2D::ResetHistory(void)
{
  /* reset if plastic */
  cout << "ResetHistory called" << endl;
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void tevp2D::Print(ostream& out) const
{
  /* inherited */
  cout << "Print called" << endl;
  FDStructMatT::Print(out);
  IsotropicT::Print(out);
  Material2DT::Print(out);
}

void tevp2D::PrintName(ostream& out) const
{
  /* inherited */
  cout << "PrintName called" << endl;
  FDStructMatT::PrintName(out);
  out << "    Thermo-Elasto-Viscoplastic\n";
}

/* spatial description */
/* modulus */
const dMatrixT& tevp2D::C_IJKL(void)
{
  /* implement modulus here */
  cout << "C_IJKL called" << endl;
  return fModulus;   // Dummy - spatial description not desired
}

/* stress */
const dSymMatrixT& tevp2D::S_IJ(void)
{
  /* implement stress here */
  cout << "S_IJ called" << endl;
  return fStress;    // Dummy - PK2 stress not desired
}

/* material description */
/* modulus */
const dMatrixT& tevp2D::c_ijkl(void)
{
  /* this function calculates the tangent modulus */
  cout << "c_ijkl called" << endl;

  /* currently computing the elastic modulus tensor */
  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
  /////// LOAD DATA?????????
  // Temporary work space arrays and matrices
  dMatrixT ata(kVoigt), lta(kVoigt);
  dArrayT pp(kVoigt), diagU(kVoigt), smlp(kVoigt);

  fModulus = ComputeDmat();   // Gives original elastic coefficient tensor
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  
  smlp = ComputeSmlp();   // Compute unnecessarily so can just access later..
  //  for (int i = 0; i < 4; i++)
  //    cout << "Smlp = " << smlp[i] << endl;
  pp = ComputePP();
  //  for (int i = 0; i < 4; i++)
  //cout << "PP = " << pp[i] << endl;
  //cout << "PP computed" << endl;
  double sb = fInternal[kSb];
  //cout << "Sb = " << sb << endl;
  double xxii = ComputeXxii();
  //cout << "xxii = " << xxii << endl;
  double ctcon = ComputeCtcon();
  //cout << "ctcon = " << ctcon << endl;
  double ebtot = ComputeEbtot();
  //cout << "Ebtot = " << ebtot << endl;
  // Compute this only so that accessor functions can get it later....

  /* calculate the tangent modulus - based on Pierce, 1987 */
  ata.Outer(pp, pp);
  lta.Outer(pp, diagU);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      /* fModulus is the tangent modulus - same as elastic modulus if the
       * gauss point has not gone plastic yet */
      fModulus(i,j) = fModulus(i,j) - ctcon * (ata(i,j) + 3.0 * El_K * Alpha_T * Chi * Xi * sb * lta(i,j));
    }
  }

  return fModulus;
}

/* stress */
const dSymMatrixT& tevp2D::s_ij(void)
{
  
  if (fRunState == GlobalT::kFormRHS)
  {
    /* implement stress here - work with Kirchoff stress, then convert back
     * to Cauchy when necessary */
    /* Allocate all state variable space on first timestep, ie time = 0 */
    cout << "s_ij called" << endl;
    const GlobalT::StateT& run_state = ContinuumElement().FEManager().RunState();
    cout << "run state = " << run_state << '\n';	
    int current_element = CurrElementNumber();
    cout << "***************** current element number = " <<  current_element << '\n';
    int ip = CurrIP();
    cout << "CurrIP = " << ip << endl;
    ElementCardT& element = CurrentElement();
    /* load data */
    LoadData(element, ip);
    
    dMatrixT kirchoff_last(3);
    dSymMatrixT cauchy_last(2);
    kirchoff_last = ArrayToMatrix(fTempKirchoff);
    cauchy_last = ArrayToSymMatrix2D(fTempCauchy);
    cout << "Past cauchy_last.  cauchy_last = \n" << cauchy_last << endl;
    cout << "Last Kirchoff Stress = \n" << kirchoff_last << endl;
    cout << "Last EffectiveStress = " << fInternal[kSb] << endl;
    cout << "Last EffectiveStrain = " << fInternal[kEb] << endl;
    cout << "Last Temperature = " << fInternal[kTemp] << endl;
    int criticalstrain = CheckCriticalStrain(element, ip);
    int checkplastic = CheckIfPlastic(element, ip);
    //cout << "Check plastic = " << checkplastic << endl;
    //cout << "CriticalStrain = " << criticalstrain << endl;
    if (criticalstrain == kFluid) {
      /* Fluid model part - if critical strain criteria is exceeded */
      ComputeGradients();   // Need determinant of deformation gradient
      double J = fFtot.Det();  // Determinant of deformation gradient
      double temp = fInternal[kTemp];   // Use the PREVIOUS temperature
      double cm = -Gamma_d * El_E * (1.0 - J + Alpha_T * (temp - Temp_0));
      cm /= (J * (1.0 - El_V));
      dMatrixT eye_cm(3);
      eye_cm = 0.0;
      eye_cm.PlusIdentity(1.0);
      eye_cm *= cm;
      fDtot *= Mu_d;
      eye_cm += fDtot;
      fStress3D = Return3DStress(eye_cm);
      dArrayT flatten_stress(kVoigt); 
      flatten_stress = MatrixToArray(fStress3D);
      fTempKirchoff = flatten_stress;
      fStress3D /= J;           // Return the Cauchy, NOT Kirchoff stress!!!
      cout << "Fluid stress computed" << endl;
    }
    else {
      /* Incremental stress update part - if critical strain criteria not
       * exceeded */
      //cout << "Still TEVP stress" << endl;
      dArrayT sig_jrate(kVoigt), dtot(kVoigt), sts_dot(kVoigt);
      sig_jrate = 0.0;
      c_ijkl();               // Need the tangent modulus
      ComputeGradients();
      double J = fFtot.Det();   // Need to convert Kirchoff to Cauchy later
      cout << "Jacobian = " << J << endl;
      /* Flatten out the Rate of Deformation tensor into Voigt notation */    
      dtot[0] = fDtot(0,0);
      dtot[1] = fDtot(1,1);
      dtot[2] = 0.0;
      dtot[3] = fDtot(0,1);   // D is symmetric - could take (1,0)
      //cout << "Rate of Deformation computed" << endl;
      //cout << "Dtot = \n" << dtot << endl;
      fModulus.Multx(dtot, sig_jrate);
      //cout << "Modulus multiplication completed" << endl;
      //cout << "sig_jrate = \n" << sig_jrate << endl;
      /* Check if plasticity has occurred yet */
      if (checkplastic == kIsPlastic)
	{
	  cout << "We're plastic!" << endl;
	  dArrayT ep_tan(kVoigt);
	  ep_tan = ComputeEP_tan();
	  sig_jrate -= ep_tan;
	}
      
      /* Add the objective part */
      //cout << "Objective part of stress starting to compute" << endl;
      sts_dot[0] = sig_jrate[0] + 2.0 * kirchoff_last(0,1) * fSpin;
      sts_dot[1] = sig_jrate[1] - 2.0 * kirchoff_last(0,1) * fSpin;
      sts_dot[2] = sig_jrate[2];
      sts_dot[3] = sig_jrate[3] - fSpin * (kirchoff_last(0,0) - kirchoff_last(1,1));
      //cout << "sts_dot = \n" << sts_dot << endl;
      double dt = GetTimeStep();
      //cout << "dt = " << dt << endl;
      kirchoff_last(0,0) += sts_dot[0] * dt;
      kirchoff_last(0,1) += sts_dot[3] * dt;
      kirchoff_last(1,0) = kirchoff_last(0,1);
      kirchoff_last(1,1) += sts_dot[1] * dt;
      kirchoff_last(2,2) += sts_dot[2] * dt;
      fStress3D = Return3DStress(kirchoff_last);
      fTempKirchoff = MatrixToArray(fStress3D);    // This is the Kirchoff stress
      
      fStress3D /= J;
    }

    fStress.ReduceFrom3D(fStress3D);     // Take only 2D stress components
    // STORE CAUCHY STRESS HERE (2D version)
    fTempCauchy = fStress;
    cout << "fStress = \n" << fStress << endl;
    cout << "fTempKirchoff = \n" << fTempKirchoff << endl;
    cout << "fTempCauchy = \n" << fTempCauchy << endl;
    /* Compute the state variables / output variables */
    fInternal[kTemp] = ComputeTemperature(element, ip);
    fInternal[kSb] = ComputeEffectiveStress();
    fInternal[kEb] = ComputeEffectiveStrain(element, ip);
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
double tevp2D::StrainEnergyDensity(void)
{
  /* compute strain energy density here */
  cout << "Computing StrainEnergyDensity" << endl;

  return fStrainEnergyDensity;
}

int tevp2D::NumOutputVariables(void) const { return kNumOutput; }
void tevp2D::OutputLabels(ArrayT<StringT>& labels) const
{
  cout << "OutputLabels called" << endl;
  /* set size */
  labels.Allocate(kNumOutput);

  /* copy labels - WHY? */
  for (int i = 0; i < kNumOutput; i++)
    labels[i] = Labels[i];
}

void tevp2D::ComputeOutput(dArrayT& output)
{
  /* Currently assuming that UpdateHistory is called before ComputeOutput */
  cout << "ComputeOutput called" << endl;
  int ip = CurrIP();
  int curr_element = CurrElementNumber();
  cout << "***************** current element number = " <<  curr_element << '\n';
  cout << "Current IP = " << ip << endl;
  ElementCardT& element = CurrentElement();
  LoadData(element, ip);
  output[0] = fInternal[kTemp];        // Temperature
  output[1] = fInternal[kEb];          // Effective strain
  output[2] = fInternal[kSb];          // Effective stress
  /* test to make sure values are right */
  cout << "fInternal[kSb] = " << fInternal[kSb] << endl;
  cout << "fInternal[kTemp] = " << fInternal[kTemp] << endl;
  cout << "fInternal[kEb] = " << fInternal[kEb] << endl;
  cout << "fTempKirchoff in ComputeOutput = \n" << fTempKirchoff << endl;
  cout << "fTempCauchy in ComputeOutput = \n" << fTempCauchy << endl;
}

/*******************************************************************
Computational routines start here
*******************************************************************/

void tevp2D::ComputeGradients(void)
{
  /* Compute the deformation gradient, rate of deformation and spin */
  /* compute deformation gradient */
  //cout << "ComputeGradients called" << endl;
  fFtot_2D = F();        
  fFtot.Rank2ExpandFrom2D(fFtot_2D);
  fFtot(2,2) = 1.0;
  //  cout << "F = \n" << fFtot << endl;
  /* compute rate of deformation */
  fF_temp.Inverse(fFtot);         // Inverse of deformation gradient  
  //cout << "fLocVel = \n" << fLocVel << endl;
  fShapes.GradU(fLocVel, fGradV_2D);    // Velocity gradient
  //  cout << "Velocity Gradient = \n" << fGradV_2D << endl;
  fGradV.Rank2ExpandFrom2D(fGradV_2D);
  fDtot.MultAB(fGradV, fF_temp, 0);  // D = dv/dx*inv(F)
  //  cout << "D = \n" << fDtot << endl;
  //  cout << "fGradV = \n" << fGradV << endl;
  /* compute spin */
  fSpin = fGradV(0,0) * fF_temp(0,1) + fGradV(0,1) * fF_temp(1,1);
  fSpin = fSpin - fGradV(1,0) * fF_temp(0,0) - fGradV(1,1) * fF_temp(1,0);
  fSpin *= .5;
  //cout << "fSpin = " << fSpin << endl;
}

dMatrixT& tevp2D::CauchyToKirchoff(dMatrixT temp_stress)
{
  /* Converts Cauchy stress to Kirchoff stress - necessary because the entire
   * stress update formulation is done in terms of Kirchoff stress */
  //cout << "CauchyToKirchoff called" << endl;
  ComputeGradients();
  double J = fFtot.Det();
  temp_stress *= J;
  fKirchoff.Swap(temp_stress);    // Copies fStress into fKirchoff, I think
  return fKirchoff;
}

double tevp2D::ComputeTemperature(ElementCardT& element, int ip)
{
  /* Compute the output temperature - 2 different methods of computing, 
   * which depends upon whether fluid model was used or not */

  /* First check to see if critical strain criteria is met */
  //cout << "ComputeTemperature called" << endl;
  double dt = GetTimeStep();
  double temp_last = fInternal[kTemp];
  int criticalstrain = CheckCriticalStrain(element, ip);
  if (criticalstrain == kFluid) {
  /* Case where fluid model was used */
    dMatrixT temp_stress(3); 
    temp_stress = ArrayToMatrix(fTempKirchoff);
    CauchyToKirchoff(temp_stress);
    ComputeGradients();
    double wpdot = fKirchoff(0,0) * fDtot(0,0) + fKirchoff(1,1) * fDtot(1,1) + 2.0 * Mu_d * pow(fDtot(0,1), 2);
    double temp_rate = Chi * Xi * wpdot;
    fTemperature = temp_rate * dt + temp_last;
    fInternal[kTemp] = fTemperature;
  }
  else {
  /* Case where fluid model was not used - viscoplasticity */
    double ebtot = GetEbtot();
    double temp_rate = Chi * Xi * ebtot;
    fTemperature = temp_rate * dt + temp_last;
    fInternal[kTemp] = fTemperature;
  }
  cout << "Computed Temperature = " << fTemperature << endl;
  return fTemperature;
}

double tevp2D::ComputeEffectiveStrain(ElementCardT& element, int ip)
{
  /* Computes the effective strain - 2 different methods of computing,
   * which depends upon whether fluid model was used or not */

  /* First check to see if critical strain criteria is met */
  //cout << "ComputeEffectiveStrain called" << endl;
  double eb_last = fInternal[kEb];
  double dt = GetTimeStep();
  int criticalstrain = CheckCriticalStrain(element, ip);

  if (criticalstrain == kFluid) {
  /* If fluid model is used (ie shear band has formed) */
    ComputeGradients();
    double temp1 = pow(fDtot(0,0), 2);
    double temp2 = pow(fDtot(1,1), 2);
    double temp3 = pow(fDtot(0,1), 2);
    double ebar = 2.0 * (temp1 + temp2) / 3.0 + 2.0 * temp3;  
    ebar = sqrt(ebar);
    fEb = ebar * dt + eb_last;
    fInternal[kEb] = fEb;
  }
  else {
  /* If fluid model / critical strain is not used */
    double ecc = ComputeEcc();

    /* access necessary data */
    double ebtot = GetEbtot();
    double ctcon = GetCtcon();
    double xxii = GetXxii();
    double ebtot_c = ebtot / (1.0 + xxii) + ctcon * ecc;
    fEb = eb_last + dt * ebtot_c;
    fInternal[kEb] = fEb;
  }
  cout << "Computed EffectiveStrain = " << fEb << endl;
  return fEb;
}

double tevp2D::ComputeEffectiveStress(void)
{
  /* Computes the effective stress - whether damage/fluid model is used
   * is apparently not relevant - it's computed the same way */
  //cout << "ComputeEffectiveStress called" << endl;
  dMatrixT temp_stress(3);
  temp_stress = ArrayToMatrix(fTempKirchoff);
  CauchyToKirchoff(temp_stress);  // Convert Cauchy stress to Kirchoff stress
  double trace = fKirchoff.Trace() / 3.0;
  double temp1 = pow(fKirchoff(0,0) - trace, 2);
  double temp2 = pow(fKirchoff(1,1) - trace, 2);
  double temp3 = pow(fKirchoff(2,2) - trace, 2);
  fSb = 1.5 * (temp1 + temp2 + temp3) + 3.0 * pow(fKirchoff(0,1), 2);
  fSb = sqrt(fSb);
  fInternal[kSb] = fSb;
  cout << "Computed EffectiveStress = " << fSb << endl;
  return fSb;
}

int tevp2D::CheckCriticalStrain(ElementCardT& element, int ip)
{
  /* Returns an indicator to determine whether critical strain criteria
   * has been met, and switch to fluid model happens next time step */
  //cout << "CheckCriticalStrain called" << endl;
  int totalIP = fNumIP;
  iArrayT& flags = element.IntegerData();
  double eb = fInternal[kEb];
  /* if already fluid, no need to check criterion */
  if (flags[ip + totalIP] == kFluid)
  {
    //cout << "Is there a problem?" << endl;
    return 1;
  }
  double ebtot = GetEbtot();  
  double ebar_cr = Epsilon_1 + (Epsilon_2 - Epsilon_1) * Epsilon_rate;
  ebar_cr /= (Epsilon_rate + ebtot);
  if (eb >= ebar_cr)
  {
    flags[ip + totalIP] = kFluid;
    fCriticalStrain = 1;        // Indicator to switch to fluid model
  }
  else
  {
    flags[ip + totalIP] = kTevp; 
    fCriticalStrain = 0;
  }
  return fCriticalStrain;
}

double tevp2D::ComputeEbtot(void)
{
  /* Computes the incremental effective strain */
  //cout << "ComputeEbtot called" << endl;
  double eb = fInternal[kEb];
  double sb = fInternal[kSb];
  double temp = fInternal[kTemp];
  if (sb <= kYieldTol)
    fEbtot = 0.0;
  else
  {
    double gsoft = Sb0 * pow(1.0 + eb/Eb0, BigN);
    gsoft *= (1.0 - Delta * (exp((temp - Temp_0)/Kappa) - 1.0));
    double reg = sb / gsoft;
    fEbtot = Eb0tot * pow(reg, Smm);
    cout << "gsoft = " << gsoft << endl;
    cout << "fEbtot = " << fEbtot << endl;
  }
  return fEbtot;
}

double tevp2D::ComputeXxii(void)
{
  /* compute Xxii - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  //cout << "ComputeXxii called" << endl;
  dArrayT pp, smlp;
  pp = GetPP();
  smlp = GetSmlp();
  double pCp2 = dArrayT::Dot(smlp, pp);
  pCp2 += Pcp;
  double sb = fInternal[kSb];
  double eb = fInternal[kEb];
  double temp = fInternal[kTemp];
  if (sb <= kYieldTol)
      fXxii = 0.0;
  else
  {
    double gsoft = Sb0 * pow(1.0 + eb/Eb0, BigN);
    gsoft *= (1.0 - Delta * (exp((temp - Temp_0)/Kappa) - 1.0));
    double reg = sb / gsoft;
    double ebtot = Eb0tot * pow(reg, Smm);
    double pE_ptau = Smm * ebtot / sb;
    double dG_Eb = BigN * gsoft / (eb + Eb0);
    double dG_T = -(Sb0 * Delta / Kappa) * (pow(1.0 + eb / Eb0, BigN) * exp((temp - Temp_0) / Kappa));
    double pE_evt = -dG_Eb * (sb / gsoft);
    double pE_Tvt = -dG_T * (sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * Chi * Xi * sb;
  
    /* Obtain the timestep */
    double dt = GetTimeStep();

    fXxii = Theta * dt * hh * pE_ptau;
  }
  cout << "fXxii = " << fXxii << endl;
  return fXxii;
}

double tevp2D::ComputeCtcon(void)
{
  /* compute Ctcon - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  //cout << "ComputeCtcon called" << endl;
  double sb = fInternal[kSb];
  double eb = fInternal[kEb];
  double temp = fInternal[kTemp];
			      
  if (sb <= kYieldTol)
    fCtcon = 0.0;
  else
  {
    dArrayT pp, smlp;
    smlp = GetSmlp();
    pp = GetPP();
    double pCp2 = dArrayT::Dot(smlp, pp);
    pCp2 += Pcp;
    double gsoft = Sb0 * pow(1.0 + eb/Eb0, BigN);
    gsoft *= (1.0 - Delta * (exp((temp - Temp_0)/Kappa) - 1.0));
    double reg = sb / gsoft;
    double ebtot = Eb0tot * pow(reg, Smm);
    double pE_ptau = Smm * ebtot / sb;
    double dG_Eb = BigN * gsoft / (eb + Eb0);
    double dG_T = -(Sb0 * Delta / Kappa) * (pow(1.0 + eb / Eb0, BigN) * exp((temp - Temp_0) / Kappa));
    double pE_evt = -dG_Eb * (sb / gsoft);
    double pE_Tvt = -dG_T * (sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * Chi * Xi * sb;
    
    /* Obtain the timestep */
    double dt = GetTimeStep();

    double xxii = Theta * dt * hh * pE_ptau;
    fCtcon = xxii / ((1.0 + xxii) * hh);
  }
  cout << "fCtcon = " << fCtcon << endl;
  return fCtcon;
}

dArrayT& tevp2D::ComputeSmlp(void)
{
  /* used in ComputePP */
  /* compute the deviatoric Kirchoff stress */
  //cout << "ComputeSmlp called" << endl;
  dMatrixT temp_stress(3); 
  temp_stress = ArrayToMatrix(fTempKirchoff);
  CauchyToKirchoff(temp_stress);
  double trace_KH = fKirchoff.Trace() / 3.0;
  fSmlp[0] = fKirchoff(0,0) - trace_KH;
  fSmlp[1] = fKirchoff(1,1) - trace_KH;
  fSmlp[2] = fKirchoff(2,2) - trace_KH;
  fSmlp[3] = fKirchoff(0,1);   // Stored in Voigt notation
  double sb = fInternal[kSb];
  if (sb <= kYieldTol)
    int blah = 0;
  else
    fSmlp *= 1.5 / sb;
  return fSmlp;
}

dArrayT& tevp2D::ComputePP(void)
{
  //cout << "ComputePP called" << endl;
  fPP = 0.0;
  dArrayT smlp(kVoigt);
  dMatrixT dmat(kVoigt);
  smlp = GetSmlp();
  dmat = GetDmat();
  dmat.Multx(smlp, fPP);
  return fPP;
}

double tevp2D::ComputeEcc(void)
{
  /* Access ecc */
  //cout << "ComputeEcc called" << endl;
  dArrayT pp(kVoigt);
  pp = GetPP();
  ComputeGradients();
  if (fInternal[kSb] <= kYieldTol)   // If hasn't yielded yet...
    fEcc = 0.0;
  else
  { 
    for (int i = 0; i < 3; i++) 
      fEcc += pp[i] * fDtot(i,i);

    fEcc += pp[3] * fDtot(0,1);
  }
  return fEcc;
}

dMatrixT& tevp2D::ComputeDmat(void)
{
  /* computes the original elastic coefficient tensor */
  //cout << "ComputeDmat called" << endl;
  fDmat = 0.0;
  double ed = El_E * (1.0 - El_V) / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double es = El_E * El_V / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double g0 = El_E / (1.0 + El_V);
  fDmat(0,0) = fDmat(1,1) = fDmat(2,2) = ed;
  fDmat(0,1) = fDmat(1,0) = fDmat(2,0) = fDmat(2,1) = es;
  fDmat(0,2) = fDmat(0,3) = fDmat(1,2) = fDmat(1,3) = 0.0;
  fDmat(2,3) = fDmat(3,0) = fDmat(3,1) = fDmat(3,2) = 0.0;
  fDmat(3,3) = g0;

  return fDmat;
}

dArrayT& tevp2D::ComputeEP_tan(void)
{
  /* computes the modulus correction if plasticity has occurred */
  cout << "ComputeEP_tan called" << endl;
  dArrayT pp(kVoigt), diagU(kVoigt);
  double ebtot = GetEbtot();
  double xxii = GetXxii();
  pp = GetPP();
  double sb = fInternal[kSb];
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  

  for (int i = 0; i < 4; i++) {
    /* EP_tan is the plastic corrector to the tangent modulus */
    fEP_tan[i] = (ebtot / (1.0 + xxii)) * (pp[i] + 3.0 * El_K * Alpha_T * Xi * Chi * sb * diagU[i]);  
  }

  return fEP_tan;
}

int tevp2D::CheckIfPlastic(ElementCardT& element, int ip)
{
  /* Checks to see if the gauss point has gone plastic yet via a
   * test on the effective stress */
  /* plastic */
  //cout << "CheckIfPlastic called" << endl;
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

void tevp2D::AllocateElement(ElementCardT& element)
{
  /* return a pointer to a new plastic element object constructed with
   * the data from element */
  /* determine storage */
  //cout << "AllocateElement called" << endl;
  int totalIP = fNumIP;                // Get total# of gauss points
  int i_size = 0;
  int d_size = 0;
  i_size += 2 * totalIP;              // 2 flags per IP:  critical strain
                                      // and check for plasticity
  d_size += kNumOutput * totalIP;     // 3 internal variables to track
  d_size += kVoigt * totalIP;         // 4 non-zero stress components:
                                      // Sig11, Sig12=Sig21, Sig22 and Sig33
  d_size += kNSD * totalIP;           // 3 2D symmetric components (Sig11, Sig12, Sig22)
  /* construct new plastic element */
  element.Allocate(i_size, d_size);

  /* first set of flags for plasticity criterion */
  for (int ip = 0; ip < totalIP; ip++)
    (element.IntegerData())[ip] = kIsElastic;
  
  /* second set of flags for critical strain / model switch criterion */
  for (int ip = totalIP; ip < 2 * totalIP; ip++)
    (element.IntegerData())[ip] = kTevp;

  element.DoubleData() = 0.0;
}

void tevp2D::LoadData(const ElementCardT& element, int ip)
{
  /* load element data for the specified integration point */
  /* check */
  //cout << "LoadData called" << endl;
  int totalIP = fNumIP;
  if (!element.IsAllocated()) throw eGeneralFail;

  int dex = ip * kVoigt;     // 4 non-zero stress components (11, 12, 22, 33)
  int dex2 = ip * kNSD;      // 3 non-zero 2D stress components (11, 12=21, 22)
  int offset = totalIP * 4;
  int offset2 = kNSD * offset / kVoigt;
  /* fetch arrays */
  dArrayT& d_array = element.DoubleData();
  fTempKirchoff.Set(kVoigt, &d_array[dex]);
  fTempCauchy.Set(kNSD, &d_array[offset + dex2]);
  fInternal.Set(kNumOutput, &d_array[offset + offset2 + ip * kNumOutput]); 
}

void tevp2D::Update(ElementCardT& element)
{
  /* get flags */
  cout << "Update called" << endl;
  iArrayT& flags = element.IntegerData();
  int totalIP = fNumIP;
  /* check if reset state (is same for all ip) */
  if (flags[0] == kReset)
  {
      flags = kIsElastic;          // don't update again
      return;
  }

  /* update plastic variables */
  for (int ip = 0; ip < totalIP; ip++)
  {
    /* fetch element data */
    LoadData(element, ip);

    /* Could make ComputeTemperature/Stress/Strain only for previous
     * timestep, then do the update here, ie fInternal += fValueComputed.
     * Then need to make the Compute functions purely for the current
     * timestep..........*/

    /* Temperature */
    //fInternal[kTemp] = ComputeTemperature(element, ip);
    
    /* Effective Stress */
    //fInternal[kSb] = ComputeEffectiveStress();
    
    /* Effective Strain */
    //fInternal[kEb] = ComputeEffectiveStrain(element, ip);
  }
}

void tevp2D::Reset(ElementCardT& element)
{
  /* resets to the last converged solution */
  /* flag to not update again */
  cout << "Reset called" << endl;
  (element.IntegerData()) = kReset;
}

dArrayT& tevp2D::MatrixToArray(dSymMatrixT StressMatrix)
{
  /* Flattens Kirchoff stress matrix into array form for internal variable
   * storage */
  //cout << "MatrixToArray called" << endl;
  fStressArray[0] = StressMatrix[0];        // Sigma 11
  fStressArray[3] = StressMatrix[5];        // Sigma 12
  fStressArray[1] = StressMatrix[1];        // Sigma 22
  fStressArray[2] = StressMatrix[2];        // Sigma 33
  
  return fStressArray;
}

dMatrixT& tevp2D::ArrayToMatrix(dArrayT StressArray)
{
  /* Expands internal variable stress array to matrix form */
  //cout << "ArrayToMatrix called" << endl;
  fStressMatrix = 0.0;
  fStressMatrix(0,0) = StressArray[0];       // Sigma 11
  fStressMatrix(1,0) = fStressMatrix(0,1) = StressArray[3];   // Sigma 12
  fStressMatrix(1,1) = StressArray[1];       // Sigma 22
  fStressMatrix(2,2) = StressArray[2];       // Sigma 33
  return fStressMatrix;
}

dSymMatrixT& tevp2D::Return3DStress(dMatrixT StressMatrix)
{
  /* Takes 3D matrix and converts to 3D symmetric matrix - necessary
   * because canned functions depend on fNumSD to convert */
  //cout << "Return3DStress called" << endl;
  fStill3D[0] = StressMatrix(0,0);      // Sigma 11
  fStill3D[1] = StressMatrix(1,1);      // Sigma 22
  fStill3D[2] = StressMatrix(2,2);      // Sigma 33
  fStill3D[3] = StressMatrix(1,2);      // Sigma 23
  fStill3D[4] = StressMatrix(0,2);      // Sigma 13
  fStill3D[5] = StressMatrix(0,1);      // Sigma 12
  return fStill3D;
}

dSymMatrixT& tevp2D::ArrayToSymMatrix2D(dArrayT StressArray)
{
  /* Because dSymMatrixT is derived from dArrayT, cannot assign fStress = StressArray */
  fSymStress2D[0] = StressArray[0];     // Sigma 11
  fSymStress2D[1] = StressArray[1];     // Sigma 22
  fSymStress2D[2] = StressArray[2];     // Sigma 33

  return fSymStress2D;
}



