/* $Id: tevp2D.cpp,v 1.8 2001-05-16 03:29:48 hspark Exp $ */
/* Implementation file for thermo-elasto-viscoplastic material subroutine */
/* Created:  Harold Park (04/04/2001) */
/* Last Updated:  Harold Park (05/10/2001) */

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
  //cout << "CONSTRUCTOR INSTANTIATED AGAIN!!" << endl;
}

/* allocate element storage */
void tevp2D::PointInitialize(void)
{
  //cout << "PointInitialize called" << endl;
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
  //cout << "UpdateHistory called" << endl;
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
  //cout << "c_ijkl called" << endl;

  /* currently computing the elastic modulus tensor */
  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
  /////// LOAD DATA?????????
  LoadData(element, ip);
  // Temporary work space arrays and matrices
  dMatrixT ata(kVoigt), lta(kVoigt);
  dArrayT diagU(kVoigt);

  ComputeDmat();   // Gives original elastic coefficient tensor
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  
  ComputeSmlp();   // Compute unnecessarily so can just access later..
  //for (int i = 0; i < 4; i++)
  //  cout << "Smlp = " << smlp[i] << endl;
  ComputePP();
  //for (int i = 0; i < 4; i++)
  //cout << "PP = " << pp[i] << endl;
  //cout << "Sb = " << sb << endl;
  //double xxii = ComputeXxii();
  //cout << "xxii = " << xxii << endl;
  //double ctcon = ComputeCtcon();
  //cout << "ctcon = " << ctcon << endl;
  //double ebtot = ComputeEbtot();
  ComputeEbtotXxiiCtcon();
  //cout << "Ebtot = " << ebtot << endl;
  // Compute this only so that accessor functions can get it later....

  /* calculate the tangent modulus - based on Pierce, 1987 */
  ata.Outer(fPP, fPP);
  lta.Outer(fPP, diagU);

  double k = fCtcon * 3.0 * El_K * Alpha_T * Chi * Xi * fInternal[kSb];
  fDmat.AddCombination(-fCtcon, ata, -k, lta);
  fModulus = fDmat;
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
    //cout << "s_ij called" << endl;
    //cout << "run state = " << run_state << '\n';	
    //cout << "***************** current element number = " <<  current_element << '\n';
    int ip = CurrIP();
    //cout << "CurrIP = " << ip << endl;
    ElementCardT& element = CurrentElement();
    /* load data */
    LoadData(element, ip);
    
    dMatrixT kirchoff_last(3);
    dSymMatrixT cauchy_last(2);
    kirchoff_last = ArrayToMatrix(fTempKirchoff);
    cauchy_last = ArrayToSymMatrix2D(fTempCauchy);
    //cout << "Past cauchy_last.  cauchy_last = \n" << cauchy_last << endl;
    //cout << "Last Kirchoff Stress = \n" << kirchoff_last << endl;
    //cout << "Last EffectiveStress = " << fInternal[kSb] << endl;
    //cout << "Last EffectiveStrain = " << fInternal[kEb] << endl;
    //cout << "Last Temperature = " << fInternal[kTemp] << endl;
    int criticalstrain = CheckCriticalStrain(element, ip);
    int checkplastic = CheckIfPlastic(element, ip);
    //cout << "Check plastic = " << checkplastic << endl;
    //cout << "CriticalStrain = " << criticalstrain << endl;
    ComputeGradients();

    if (criticalstrain == kFluid) {
      /* Fluid model part - if critical strain criteria is exceeded */
      //cout << "FAILURE MODEL SWITCH!!" << endl;
      //ComputeGradients();   // Need determinant of deformation gradient
      double J = fFtot.Det();  // Determinant of deformation gradient
      double temp = fInternal[kTemp];   // Use the PREVIOUS temperature
      double cm = -Gamma_d * El_E * (1.0 - J + Alpha_T * (temp - Temp_0));
      cm /= (J * (1.0 - El_V));
      dMatrixT eye_cm(3), dtemp(3);
      eye_cm = dtemp = 0.0;
      eye_cm.PlusIdentity(1.0);
      eye_cm *= cm;
      dtemp = fDtot;
      dtemp *= Mu_d;
      eye_cm += dtemp;
      fStress3D = Return3DStress(eye_cm);
      fTempKirchoff = MatrixToArray(fStress3D);
      fStress3D /= J;           // Return the Cauchy, NOT Kirchoff stress!!!
    }  
    else {
      /* Incremental stress update part - if critical strain criteria not
       * exceeded */
      //cout << "Still TEVP stress" << endl;
      dArrayT sig_jrate(kVoigt), dtot(kVoigt), sts_dot(kVoigt);
      sig_jrate = 0.0;
      c_ijkl();               // Need the tangent modulus
      //ComputeGradients();
      double J = fFtot.Det();   // Need to convert Kirchoff to Cauchy later
      //cout << "Jacobian = " << J << endl;
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
	sig_jrate -= ComputeEP_tan();
      
      /* Add the objective part */
      //cout << "Objective part of stress starting to compute" << endl;
      sts_dot[0] = sig_jrate[0] + 2.0 * kirchoff_last(0,1) * fSpin;
      sts_dot[1] = sig_jrate[1] - 2.0 * kirchoff_last(0,1) * fSpin;
      sts_dot[2] = sig_jrate[2];
      sts_dot[3] = sig_jrate[3] - fSpin * (kirchoff_last(0,0) - kirchoff_last(1,1));
      //cout << "sts_dot = \n" << sts_dot << endl;
      //cout << "dt = " << dt << endl;
      kirchoff_last(0,0) += sts_dot[0] * fDt;
      kirchoff_last(0,1) += sts_dot[3] * fDt;
      kirchoff_last(1,0) = kirchoff_last(0,1);
      kirchoff_last(1,1) += sts_dot[1] * fDt;
      kirchoff_last(2,2) += sts_dot[2] * fDt;
      fStress3D = Return3DStress(kirchoff_last);
      fTempKirchoff = MatrixToArray(fStress3D);    // This is the Kirchoff stress
      
      fStress3D /= J;
    }

    fStress.ReduceFrom3D(fStress3D);     // Take only 2D stress components
    // STORE CAUCHY STRESS HERE (2D version)
    fTempCauchy = fStress;
    //cout << "fStress = \n" << fStress << endl;
    //cout << "fTempKirchoff = \n" << fTempKirchoff << endl;
    //cout << "fTempCauchy = \n" << fTempCauchy << endl;
    /* Compute the state variables / output variables */
    fInternal[kSb] = ComputeEffectiveStress();
    fInternal[kTemp] = ComputeTemperature();
    fInternal[kEb] = ComputeEffectiveStrain();
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
  //cout << "ComputeOutput called" << endl;
  int ip = CurrIP();
  //cout << "Current IP = " << ip << endl;
  ElementCardT& element = CurrentElement();
  //int criticalstrain = CheckCriticalStrain(element, ip);
  //int checkplastic = CheckIfPlastic(element, ip);
  //ComputeGradients();
  //double ebtot = GetEbtot();
  //double xxii = GetXxii();
  //dArrayT smlp = GetSmlp();
  //dArrayT pp = GetPP();
  //double ctcon = GetCtcon();
  LoadData(element, ip);
  output[0] = fInternal[kTemp];        // Temperature
  output[1] = fInternal[kEb];          // Effective strain
  output[2] = fInternal[kSb];          // Effective stress
  /* test to make sure values are right */
  //cout << "CriticalStrain = " << criticalstrain << endl;
  //cout << "CheckPlastic = " << checkplastic << endl;
  //cout << "Ebtot = " << ebtot << endl;
  //cout << "PP = " << pp << endl;
  //cout << "smlp = " << smlp << endl;
  //cout << "ctcon = " << ctcon << endl;
  //cout << "xxii = " << xxii << endl;
  //cout << "fInternal[kSb] = " << fInternal[kSb] << endl;
  //cout << "fInternal[kTemp] = " << fInternal[kTemp] << endl;
  //cout << "fInternal[kEb] = " << fInternal[kEb] << endl;
  //cout << "Kirchoff Stress in ComputeOutput = \n" << fTempKirchoff << '\n';
  //cout << "fTempCauchy in ComputeOutput = \n" << fTempCauchy << endl;
  //cout << "Rate of Deformation = \n" << fDtot << '\n';
  //cout << "Deformation Gradient = \n" << fFtot << '\n';
  //cout << "Velocity Strain = \n" << fGradV << '\n';
  //cout << "Inverse of Deformation Gradient = \n" << fF_temp << '\n';
}

/*******************************************************************
Computational routines start here
*******************************************************************/

void tevp2D::ComputeGradients(void)
{
  /* Compute the deformation gradient, rate of deformation and spin */
  /* compute deformation gradient */
  //cout << "ComputeGradients called" << endl;
  double& tempspin = fSpin;
  dMatrixT& tempf = fFtot;
  dMatrixT& tempd = fDtot;
  tempf = tempd = 0.0;

  fFtot_2D = F();        
  tempf.Rank2ExpandFrom2D(fFtot_2D);
  tempf(2,2) = 1.0;
  //cout << "F = \n" << fFtot << endl;
  /* compute rate of deformation */
  fF_temp.Inverse(tempf);         // Inverse of deformation gradient  
  //cout << "fLocVel = \n" << fLocVel << '\n';
  //cout << "Inverse of deformation gradient = \n" << fF_temp << '\n';
  fShapes.GradU(fLocVel, fGradV_2D);    // Velocity gradient
  //cout << "Velocity Gradient = \n" << fGradV_2D << '\n';
  fGradV.Rank2ExpandFrom2D(fGradV_2D);
  //fDtot.MultAB(fGradV, fF_temp, 0);  // D = dv/dx*inv(F)
  tempd(0,0) = fGradV(0,0) * fF_temp(0,0) + fGradV(0,1) * fF_temp(1,0);
  tempd(1,1) = fGradV(1,0) * fF_temp(0,1) + fGradV(1,1) * fF_temp(1,1);
  tempd(2,2) = tempd(1,2) = tempd(2,1) = 0.0;
  tempd(0,1) = tempd(1,0) = fGradV(0,0) * fF_temp(0,1) + fGradV(0,1) * fF_temp(1,1) + fGradV(1,0) * fF_temp(0,0) + fGradV(1,1) * fF_temp(1,0);
  tempd(0,1) *= .5;
  tempd(1,0) *= .5;
  //cout << "D = \n" << fDtot << '\n';
  //cout << "fGradV = \n" << fGradV << '\n';
  /* compute spin */
  tempspin = fGradV(0,0) * fF_temp(0,1) + fGradV(0,1) * fF_temp(1,1);
  tempspin = tempspin - fGradV(1,0) * fF_temp(0,0) - fGradV(1,1) * fF_temp(1,0);
  tempspin *= .5;
  //cout << "fSpin = " << fSpin << '\n';
}

double tevp2D::ComputeTemperature(void)
{
  /* Compute the output temperature - 2 different methods of computing, 
   * which depends upon whether fluid model was used or not */

  /* First check to see if critical strain criteria is met */
  //cout << "ComputeTemperature called" << endl;
  double temp_last = fInternal[kTemp];

  if (fCriticalStrain == kFluid) {
  /* Case where fluid model was used */
    dMatrixT temp_stress(3); 
    temp_stress = ArrayToMatrix(fTempKirchoff);
    //ComputeGradients();
    double wpdot = temp_stress(0,0) * fDtot(0,0) + temp_stress(1,1) * fDtot(1,1) + 2.0 * Mu_d * pow(fDtot(0,1), 2);
    double temp_rate = Chi * Xi * wpdot;
    fTemperature = temp_rate * fDt + temp_last;
    if (temp_rate < 0.0)
      cout << "NEGATIVE TEMPERATURE RATE:  FLUID" << '\n';
    fInternal[kTemp] = fTemperature;
  }
  else {
  /* Case where fluid model was not used - viscoplasticity */
    double temp_rate = Chi * Xi * fEbtot * fInternal[kSb];
    fTemperature = temp_rate * fDt + temp_last;
    if (temp_rate < 0.0)
      cout << "NEGATIVE TEMPERATURE RATE:  VISCO" << '\n';
    fInternal[kTemp] = fTemperature;
  }
  //cout << "Computed Temperature = " << fTemperature << endl;
  if (fTemperature < 293.0)
    cout << "TEMPERATURE < 293K!!!" << '\n';
  return fTemperature;
}

double tevp2D::ComputeEffectiveStrain(void)
{
  /* Computes the effective strain - 2 different methods of computing,
   * which depends upon whether fluid model was used or not */

  /* First check to see if critical strain criteria is met */
  //cout << "ComputeEffectiveStrain called" << endl;
  double eb_last = fInternal[kEb];

  if (fCriticalStrain == kFluid) {
  /* If fluid model is used (ie shear band has formed) */
    //ComputeGradients();
    double temp1 = pow(fDtot(0,0), 2);
    double temp2 = pow(fDtot(1,1), 2);
    double temp3 = pow(fDtot(0,1), 2);
    double ebar = 2.0 * (temp1 + temp2) / 3.0 + 2.0 * temp3;  
    ebar = sqrt(ebar);
    fEb = ebar * fDt + eb_last;
    if (ebar < 0.0)
      cout << "NEGATIVE EFFECTIVE STRAIN RATE:  FLUID" << '\n';
    fInternal[kEb] = fEb;
  }
  else {
  /* If fluid model / critical strain is not used */
    double ecc = ComputeEcc();

    /* access necessary data */
    double ebtot_c = fEbtot / (1.0 + fXxii) + fCtcon * ecc;
    fEb = eb_last + fDt * ebtot_c;
    //if (ebtot_c < 0.0)
    //  cout << "NEGATIVE EFFECTIVE STRAIN RATE:  VISCO" << '\n';
    fInternal[kEb] = fEb;
  }
  //cout << "Computed EffectiveStrain = " << fEb << endl;
  if (fEb < 0.0)
    cout << "TOTAL NEGATIVE EFFECTIVE STRAIN" << '\n';
  return fEb;
}

double tevp2D::ComputeEffectiveStress(void)
{
  /* Computes the effective stress - whether damage/fluid model is used
   * is apparently not relevant - it's computed the same way */
  //cout << "ComputeEffectiveStress called" << endl;
  dMatrixT temp_stress(3);
  temp_stress = ArrayToMatrix(fTempKirchoff);
  double trace = temp_stress.Trace() / 3.0;
  double temp1 = pow(temp_stress(0,0) - trace, 2);
  double temp2 = pow(temp_stress(1,1) - trace, 2);
  double temp3 = pow(temp_stress(2,2) - trace, 2);
  fSb = 1.5 * (temp1 + temp2 + temp3) + 3.0 * pow(temp_stress(0,1), 2);
  fSb = sqrt(fSb);
  fInternal[kSb] = fSb;
  //cout << "Computed EffectiveStress = " << fSb << '\n';
  if (fSb < 0.0)
    cout << "NEGATIVE EFFECTIVE STRESS" << '\n';
  return fSb;
}

int tevp2D::CheckCriticalStrain(const ElementCardT& element, int ip)
{
  /* Returns an indicator to determine whether critical strain criteria
   * has been met, and switch to fluid model happens next time step */
  //cout << "CheckCriticalStrain called" << endl;
  iArrayT& flags = element.IntegerData();
  int& criticalstrain = fCriticalStrain;
  /* if already fluid, no need to check criterion */
  if (flags[ip + fNumIP] == kFluid)
  {
    criticalstrain = 1;
    return criticalstrain;
  }
  double ebar_cr = Epsilon_1 + (Epsilon_2 - Epsilon_1) * Epsilon_rate;
  ebar_cr /= (Epsilon_rate + fEbtot);
  if (fInternal[kEb] >= ebar_cr)
  {
    flags[ip + fNumIP] = kFluid;
    criticalstrain = 1;        // Indicator to switch to fluid model
  }
  else
  {
    flags[ip + fNumIP] = kTevp; 
    criticalstrain = 0;
  }
  return criticalstrain;
}

void tevp2D::ComputeEbtotXxiiCtcon(void)
{
  /* compute Ctcon - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  //cout << "ComputeCtcon called" << endl;
  double& ebtot = fEbtot;
  double& xxii = fXxii;
  double& ctcon = fCtcon;
  double sb = fInternal[kSb];
  double eb = fInternal[kEb];
  double temp = fInternal[kTemp];
			      
  if (sb <= kYieldTol)
  {
    ebtot = 0.0;
    xxii = 0.0;
    ctcon = 0.0;
  }
  else
  {
    double pCp2 = dArrayT::Dot(fSmlp, fPP);
    pCp2 += Pcp;
    double gsoft = Sb0 * pow(1.0 + eb/Eb0, BigN);
    gsoft *= (1.0 - Delta * (exp((temp - Temp_0)/Kappa) - 1.0));
    double reg = sb / gsoft;
    ebtot = Eb0tot * pow(reg, Smm);
    double pE_ptau = Smm * ebtot / sb;
    double dG_Eb = BigN * gsoft / (eb + Eb0);
    double dG_T = -(Sb0 * Delta / Kappa) * (pow(1.0 + eb / Eb0, BigN) * exp((temp - Temp_0) / Kappa));
    double pE_evt = -dG_Eb * (sb / gsoft);
    double pE_Tvt = -dG_T * (sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * Chi * Xi * sb;
    
    xxii = Theta * fDt * hh * pE_ptau;
    ctcon = xxii / ((1.0 + xxii) * hh);
  }
  if (ctcon < 0.0)
    cout << "NEGATIVE CTCON!" << '\n';

}

void tevp2D::ComputeSmlp(void)
{
  /* used in ComputePP */
  /* compute the deviatoric Kirchoff stress */
  //cout << "ComputeSmlp called" << endl;
  dArrayT& smlp = fSmlp;
  dMatrixT temp_stress(3); 
  temp_stress = ArrayToMatrix(fTempKirchoff);
  double trace_KH = temp_stress.Trace() / 3.0;
  smlp[0] = temp_stress(0,0) - trace_KH;
  smlp[1] = temp_stress(1,1) - trace_KH;
  smlp[2] = temp_stress(2,2) - trace_KH;
  smlp[3] = temp_stress(0,1);   // Stored in Voigt notation
  if (fInternal[kSb] <= kYieldTol) 
    smlp = 0.0;
  else
    smlp *= 1.5 / fInternal[kSb];

}

void tevp2D::ComputePP(void)
{
  //cout << "ComputePP called" << endl;
  dArrayT& pp = fPP;
  fDmat.Multx(fSmlp, pp);
}

double tevp2D::ComputeEcc(void)
{
  /* Access ecc */
  //cout << "ComputeEcc called" << endl;
  fEcc = 0.0;
  //ComputeGradients();
  if (fInternal[kSb] <= kYieldTol)   // If hasn't yielded yet...
    fEcc = 0.0;
  else
  { 
    for (int i = 0; i < 3; i++) 
      fEcc += fPP[i] * fDtot(i,i);

    fEcc += fPP[3] * fDtot(0,1);
  }
  //if (fEcc < 0.0)
  //cout << "NEGATIVE ECC" << '\n';
  return fEcc;
}

void tevp2D::ComputeDmat(void)
{
  /* computes the original elastic coefficient tensor */
  //cout << "ComputeDmat called" << endl;
  dMatrixT& dmat = fDmat;
  double ed = El_E * (1.0 - El_V) / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double es = El_E * El_V / ((1.0 + El_V) * (1.0 - 2.0 * El_V));
  double g0 = El_E / (1.0 + El_V);
  dmat(0,0) = dmat(1,1) = dmat(2,2) = ed;
  dmat(0,1) = dmat(1,0) = dmat(2,0) = dmat(2,1) = es;
  dmat(0,2) = dmat(0,3) = dmat(1,2) = dmat(1,3) = 0.0;
  dmat(2,3) = dmat(3,0) = dmat(3,1) = dmat(3,2) = 0.0;
  dmat(3,3) = g0;

}

dArrayT& tevp2D::ComputeEP_tan(void)
{
  /* computes the modulus correction if plasticity has occurred */
  //cout << "ComputeEP_tan called" << endl;
  fEP_tan = 0.0;
  dArrayT diagU(kVoigt);
  diagU[0] = diagU[1] = diagU[2] = 1.0;
  diagU[3] = 0.0;  

  for (int i = 0; i < 4; i++) {
    /* EP_tan is the plastic corrector to the tangent modulus */
    fEP_tan[i] = (fEbtot / (1.0 + fXxii)) * (fPP[i] + 3.0 * El_K * Alpha_T * Xi * Chi * fInternal[kSb] * diagU[i]);  
  }

  return fEP_tan;
}

int tevp2D::CheckIfPlastic(const ElementCardT& element, int ip)
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
  int i_size = 0;
  int d_size = 0;
  i_size += 2 * fNumIP;              // 2 flags per IP:  critical strain
                                      // and check for plasticity
  d_size += kNumOutput * fNumIP;     // 3 internal variables to track
  d_size += kVoigt * fNumIP;         // 4 non-zero stress components:
                                      // Sig11, Sig12=Sig21, Sig22 and Sig33
  d_size += kNSD * fNumIP;           // 3 2D symmetric components (Sig11, Sig12, Sig22)
  /* construct new plastic element */
  element.Allocate(i_size, d_size);

  /* first set of flags for plasticity criterion */
  for (int ip = 0; ip < fNumIP; ip++)
    (element.IntegerData())[ip] = kIsElastic;
  
  /* second set of flags for critical strain / model switch criterion */
  for (int ip = fNumIP; ip < 2 * fNumIP; ip++)
    (element.IntegerData())[ip] = kTevp;

  element.DoubleData() = 0.0;
}

void tevp2D::LoadData(const ElementCardT& element, int ip)
{
  /* load element data for the specified integration point */
  /* check */
  //cout << "LoadData called" << endl;
  if (!element.IsAllocated()) throw eGeneralFail;

  int dex = ip * kVoigt;     // 4 non-zero stress components (11, 12, 22, 33)
  int dex2 = ip * kNSD;      // 3 non-zero 2D stress components (11, 12=21, 22)
  int offset = fNumIP * 4;
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
  //cout << "Update called" << endl;
  iArrayT& flags = element.IntegerData();
  /* check if reset state (is same for all ip) */
  if (flags[0] == kReset)
  {
      flags = kIsElastic;          // don't update again
      return;
  }
}

void tevp2D::Reset(ElementCardT& element)
{
  /* resets to the last converged solution */
  /* flag to not update again */
  cout << "Reset called" << endl;
  (element.IntegerData()) = kReset;
}

dArrayT& tevp2D::MatrixToArray(const dSymMatrixT& StressMatrix)
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

dMatrixT& tevp2D::ArrayToMatrix(const dArrayT& StressArray)
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

dSymMatrixT& tevp2D::Return3DStress(const dMatrixT& StressMatrix)
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

dSymMatrixT& tevp2D::ArrayToSymMatrix2D(const dArrayT& StressArray)
{
  /* Because dSymMatrixT is derived from dArrayT, cannot assign fStress = StressArray */
  fSymStress2D[0] = StressArray[0];     // Sigma 11
  fSymStress2D[1] = StressArray[1];     // Sigma 22
  fSymStress2D[2] = StressArray[2];     // Sigma 33

  return fSymStress2D;
}



