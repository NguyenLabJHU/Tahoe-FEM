/* Implementation file for thermo-elasto-viscoplastic material subroutine */
/* Created:  Harold Park (04/04/2001) */
/* Last Updated:  Harold Park (04/30/2001) */

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
  fStress3D(3),
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
  fStressArray(kVoigt),
  fSmlp(kVoigt)
  
{
  /* initialize material constants */
  el_E = 2.0E11;
  el_V = .30;
  el_K = el_E / (3.0 * (1.0 - 2.0 * el_V));
  el_G = el_E / (2.0 * (1.0 + el_V));
  Sb0 = 2.0E9;
  rho0 = 7830.0;
  Eb0 = Sb0 / el_E;
  Eb0tot = .001;
  bigN = .01;                // strain hardening exponent
  smm = 70.0;                // rate sensitivity parameter

  /* initialize temperature parameters */
  Temp_0 = 293.0;
  alpha_T = 11.2E-6;
  delta = .8;
  theta = .5;
  kappa = 500.0;  
  Cp = 448.0;
  chi = .9;
  ecc = 0.0;
  pCp = 3.0 * el_G;

  /* initialize damage parameters */
  epsilon_1 = 4.0 * Eb0;
  epsilon_2 = .3;
  epsilon_rate = 4.0E4;
  gamma_d = .002;
  mu_d = 500.0;
  sigCr = 6.0 * Sb0;

  /* used in temperature update */
  xi = 1.0 / (rho0 * Cp);
  cout << "CONSTRUCTOR INSTANTIATED AGAIN!!" << endl;
}

/* allocate element storage */
void tevp2D::PointInitialize(void)
{
	/* first ip only */
	if (CurrIP() == 0) AllocateElement(CurrentElement());
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
  LoadData(element, ip);
  // Temporary work space arrays and matrices
  dMatrixT AtA(kVoigt), LtA(kVoigt);
  dArrayT PP(kVoigt), diagU(kVoigt), Smlp(kVoigt);

  fModulus = ComputeDmat();   // Gives original elastic coefficient tensor
  diagU[0] = 1.0;
  diagU[1] = 1.0;
  diagU[2] = 1.0;
  diagU[3] = 0.0;  
  Smlp = ComputeSmlp();   // Compute unnecessarily so can just access later..
  //  for (int i = 0; i < 4; i++)
  //    cout << "Smlp = " << Smlp[i] << endl;
  PP = ComputePP();
  //  for (int i = 0; i < 4; i++)
  //cout << "PP = " << PP[i] << endl;
  //cout << "PP computed" << endl;
  double Sb = fInternal[kSb];
  //cout << "Sb = " << Sb << endl;
  double xxii = ComputeXxii();
  //cout << "xxii = " << xxii << endl;
  double ctcon = ComputeCtcon();
  //cout << "ctcon = " << ctcon << endl;
  double Ebtot = ComputeEbtot();
  //cout << "Ebtot = " << Ebtot << endl;
  // Compute this only so that accessor functions can get it later....

  /* calculate the tangent modulus - based on Pierce, 1987 */
  AtA.Outer(PP, PP);
  LtA.Outer(PP, diagU);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      /* fModulus is the tangent modulus - same as elastic modulus if the
       * gauss point has not gone plastic yet */
      fModulus(i,j) = fModulus(i,j) - ctcon * (AtA(i,j) + 3.0 * el_K * alpha_T * chi * xi * Sb * LtA(i,j));
    }
  }

  return fModulus;
}

/* stress */
const dSymMatrixT& tevp2D::s_ij(void)
{
  /* implement stress here - work with Kirchoff stress, then convert back
   * to Cauchy when necessary */
  /* Allocate all state variable space on first timestep, ie time = 0 */
  cout << "s_ij called" << endl;
  const FEManagerT& fe_man = ContinuumElement().FEManager();
  //fTime = fe_man.Time();
  //fDt = fe_man.TimeStep();
  int ip = CurrIP();
  cout << "CurrIP = " << ip << endl;
  ElementCardT& element = CurrentElement();
  //for (int i = 0; i < 2*kVoigt; i++)
  //cout << "fLocVel = " << fLocVel[i] << endl;
//  if (fTime == fDt) // first ip only
// {
    /* initialize all state variable storage space during first timestep */
//  AllocateElement(element);
//    LoadData(element, ip);
//  }
//  else
// {
    //cout << "Else hit..." << endl;
    //LoadData(element, ip);
    //cout << "Past LoadData" << endl;
// }
 
 	/* load data */
 	LoadData(element, ip);
 
  dMatrixT Stress_Last(3);
  Stress_Last = ArrayToMatrix(fTempStress);
//  for (int i = 0; i < 3; i++)
//    for (int j = 0; j < 3; j++)
//      cout << "Stress_Last = " << Stress_Last(i,j) << endl;
  cout << "Stress_Last = \n" << Stress_Last << endl;
  
  cout << "Last EffectiveStress = " << fInternal[kSb] << endl;
  cout << "Last EffectiveStrain = " << fInternal[kEb] << endl;
  cout << "Last Temperature = " << fInternal[kTemp] << endl;
  int CriticalStrain = CheckCriticalStrain(element, ip);
  int CheckPlastic = CheckIfPlastic(element, ip);
  //cout << "Check plastic = " << CheckPlastic << endl;
  //cout << "CriticalStrain = " << CriticalStrain << endl;
  if (CriticalStrain == kFluid) {
    /* Fluid model part - if critical strain criteria is exceeded */
    ComputeGradients();   // Need determinant of deformation gradient
    cout << "After ComputeGradients?" << endl;
    double J = fFtot.Det();  // Determinant of deformation gradient
    double Temp = fInternal[kTemp];   // Use the PREVIOUS temperature
    double cm = -gamma_d * el_E * (1.0 - J + alpha_T * (Temp - Temp_0));
    cm /= (J * (1.0 - el_V));
    dMatrixT eye_cm(3);
    eye_cm = 0.0;
    eye_cm.PlusIdentity(1.0);
    eye_cm *= cm;
    fDtot *= mu_d;
    eye_cm += fDtot;
    
    //TEMP
    cout << "\n tevp2D::s_ij: materiak went fluid and there is a dimension mismatch here\n"
         <<   "     between fStress and eye_cm" << endl;
    
    fStress.FromMatrix(eye_cm);
    dArrayT Flatten_Stress(kVoigt); 
    Flatten_Stress = MatrixToArray(fStress);
    fTempStress = Flatten_Stress;
    fStress /= J;           // Return the Cauchy, NOT Kirchoff stress!!!
    cout << "Stress computed" << endl;
  }
  else {
    /* Incremental stress update part - if critical strain criteria not
     * exceeded */
    //cout << "Still TEVP stress" << endl;
    dArrayT sig_jrate(kVoigt), Dtot(kVoigt), sts_dot(kVoigt);
    sig_jrate = 0.0;
    c_ijkl();               // Need the tangent modulus
    ComputeGradients();
    double J = fFtot.Det();   // Need to convert Kirchoff to Cauchy later
    /* Flatten out the Rate of Deformation tensor into Voigt notation */    
    Dtot[0] = fDtot(0,0);
    Dtot[1] = fDtot(1,1);
    Dtot[2] = 0.0;
    Dtot[3] = fDtot(0,1);   // D is symmetric - could take (1,0)
    //cout << "Rate of Deformation computed" << endl;
    //for (int i = 0; i < 4; i++)
    //cout << "Dtot = " << Dtot[i] << endl;
    fModulus.Multx(Dtot, sig_jrate);
    //cout << "Modulus multiplication completed" << endl;
    //for (int i = 0; i < 4; i++)
    //cout << "sig_jrate = " << sig_jrate[i] << endl;
    /* Check if plasticity has occurred yet */
    if (CheckPlastic == kIsPlastic)
    {
      cout << "We're plastic!" << endl;
      dArrayT EP_tan(kVoigt);
      EP_tan = ComputeEP_tan();
      sig_jrate -= EP_tan;
    }

    /* Add the objective part */
    //cout << "Objective part of stress starting to compute" << endl;
    sts_dot[0] = sig_jrate[0] + 2.0 * Stress_Last(0,1) * fSpin;
    sts_dot[1] = sig_jrate[1] - 2.0 * Stress_Last(0,1) * fSpin;
    sts_dot[2] = sig_jrate[2];
    sts_dot[3] = sig_jrate[3] - fSpin * (Stress_Last(0,0) - Stress_Last(1,1));
    //for (int i = 0; i < 4; i++)
    //cout << "sts_dot = " << sts_dot[i] << endl;
    double dt = GetTimeStep();
    //cout << "dt = " << dt << endl;
    Stress_Last(0,0) += sts_dot[0] * dt;
    Stress_Last(0,1) += sts_dot[3] * dt;
    Stress_Last(1,0) = Stress_Last(0,1);
    Stress_Last(1,1) += sts_dot[1] * dt;
    Stress_Last(2,2) += sts_dot[2] * dt;
    fStress3D = Return3DStress(Stress_Last);
    fTempStress = MatrixToArray(fStress3D);

    fStress3D /= J;
  }

  fStress.ReduceFrom3D(fStress3D);     // Take only 2D stress components
  //  for (int i = 0; i < 3; i++)
  //cout << "fStress = " << fStress[i] << endl;
  /* Compute the state variables / output variables */
  fInternal[kTemp] = ComputeTemperature(element, ip);
  fInternal[kSb] = ComputeEffectiveStress();
  fInternal[kEb] = ComputeEffectiveStrain(element, ip);

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
  ElementCardT& element = CurrentElement();
  LoadData(element, ip);
  output[0] = fInternal[kTemp];        // Temperature
  output[1] = fInternal[kEb];          // Effective strain
  output[2] = fInternal[kSb];          // Effective stress
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
  /* compute rate of deformation */
  fF_temp.Inverse(fFtot);         // Inverse of deformation gradient  
  //for (int i = 0; i < 8; i++)
  //cout << "fLocVel = " << fLocVel[i] << endl;
  fShapes.GradU(fLocVel, fGradV_2D);    // Velocity gradient
  //for (int i = 0; i < 2; i++)
  //for (int j = 0; j < 2; j++)
  //  cout << "Velocity Gradient = " << fGradV_2D(i,j) << endl;
  fGradV.Rank2ExpandFrom2D(fGradV_2D);
  fDtot.MultAB(fGradV, fF_temp, 0);  // D = dv/dx*inv(F)

  /* compute spin */
  fSpin = fGradV(0,0) * fF_temp(0,1) + fGradV(0,1) * fF_temp(1,1);
  fSpin = fSpin - fGradV(1,0) * fF_temp(0,0) - fGradV(1,1) * fF_temp(1,0);
  fSpin *= .5;
}

dMatrixT& tevp2D::CauchyToKirchoff(dMatrixT Temp_Stress)
{
  /* Converts Cauchy stress to Kirchoff stress - necessary because the entire
   * stress update formulation is done in terms of Kirchoff stress */
  //cout << "CauchyToKirchoff called" << endl;
  ComputeGradients();
  double J = fFtot.Det();
  Temp_Stress *= J;
  fKirchoff.Swap(Temp_Stress);    // Copies fStress into fKirchoff, I think
  return fKirchoff;
}

double tevp2D::ComputeTemperature(ElementCardT& element, int ip)
{
  /* Compute the output temperature - 2 different methods of computing, 
   * which depends upon whether fluid model was used or not */

  /* First check to see if critical strain criteria is met */
  //cout << "ComputeTemperature called" << endl;
  double dt = GetTimeStep();
  cout << "dt within ComputeTemperature = " << dt << endl;
  double Temp_last = fInternal[kTemp];
  int CriticalStrain = CheckCriticalStrain(element, ip);
  if (CriticalStrain == kFluid) {
  /* Case where fluid model was used */
    dMatrixT temp_stress(3); 
    temp_stress = ArrayToMatrix(fTempStress);
    CauchyToKirchoff(temp_stress);
    ComputeGradients();
    double Wpdot = fKirchoff(0,0) * fDtot(0,0) + fKirchoff(1,1) * fDtot(1,1) + 2.0 * mu_d * pow(fDtot(0,1), 2);
    double Temp_rate = chi * xi * Wpdot;
    fTemperature = Temp_rate * dt + Temp_last;
    fInternal[kTemp] = fTemperature;
  }
  else {
  /* Case where fluid model was not used - viscoplasticity */
    double Ebtot = GetEbtot();
    cout << "Ebtot within ComputeTemperature = " << Ebtot << endl;
    double temp_rate = chi * xi * Ebtot;
    fTemperature = temp_rate * dt + Temp_last;
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
  double Eb_last = fInternal[kEb];
  double dt = GetTimeStep();
  int CriticalStrain = CheckCriticalStrain(element, ip);

  if (CriticalStrain == kFluid) {
  /* If fluid model is used (ie shear band has formed) */
    ComputeGradients();
    double temp1 = pow(fDtot(0,0), 2);
    double temp2 = pow(fDtot(1,1), 2);
    double temp3 = pow(fDtot(0,1), 2);
    double Ebar = 2.0 * (temp1 + temp2) / 3.0 + 2.0 * temp3;  
    Ebar = sqrt(Ebar);
    fEb = Ebar * dt + Eb_last;
    fInternal[kEb] = fEb;
  }
  else {
  /* If fluid model / critical strain is not used */
    double ecc = ComputeEcc();

    /* access necessary data */
    double Ebtot = GetEbtot();
    double ctcon = GetCtcon();
    double xxii = GetXxii();
    double Ebtot_c = Ebtot / (1.0 + xxii) + ctcon * ecc;
    fEb = Eb_last + dt * Ebtot_c;
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
  temp_stress = ArrayToMatrix(fTempStress);
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
  int TotalIP = fNumIP;
  iArrayT& flags = element.IntegerData();
  /* if already fluid, no need to check criterion */
  if (flags[ip + TotalIP] == kFluid)
  {
    //cout << "Is there a problem?" << endl;
    return 1;
  }
  double Ebtot = GetEbtot();  
  double Ebar_cr = epsilon_1 + (epsilon_2 - epsilon_1) * epsilon_rate;
  Ebar_cr /= (epsilon_rate + Ebtot);
  if (fEb >= Ebar_cr)
  {
    flags[ip + TotalIP] = kFluid;
    fCriticalStrain = 1;        // Indicator to switch to fluid model
  }
  else
  {
    flags[ip + TotalIP] = kTevp; 
    fCriticalStrain = 0;
  }
  return fCriticalStrain;
}

double tevp2D::ComputeEbtot(void)
{
  /* Computes the incremental effective strain */
  //cout << "ComputeEbtot called" << endl;
  double Eb = fInternal[kEb];
  double Sb = fInternal[kSb];
  double Temp = fInternal[kTemp];
  if (Sb <= kYieldTol)
    fEbtot = 0.0;
  else
  {
    double gsoft = Sb0 * pow(1.0 + Eb/Eb0, bigN);
    gsoft *= (1.0 - delta * (exp((Temp - Temp_0)/kappa) - 1.0));
    double reg = Sb / gsoft;
    fEbtot = Eb0tot * pow(reg, smm);
  }
  return fEbtot;
}

double tevp2D::ComputeXxii(void)
{
  /* compute Xxii - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  //cout << "ComputeXxii called" << endl;
  dArrayT PP, Smlp;
  PP = GetPP();
  Smlp = GetSmlp();
  double pCp2 = dArrayT::Dot(Smlp, PP);
  pCp2 += pCp;
  double Sb = fInternal[kSb];
  double Eb = fInternal[kEb];
  double Temp = fInternal[kTemp];
  if (Sb <= kYieldTol)
      fXxii = 0.0;
  else
  {
    double gsoft = Sb0 * pow(1.0 + Eb/Eb0, bigN);
    gsoft *= (1.0 - delta * (exp((Temp - Temp_0)/kappa) - 1.0));
    double reg = Sb / gsoft;
    double Ebtot = Eb0tot * pow(reg, smm);
    double pE_ptau = smm * Ebtot / Sb;
    double dG_Eb = bigN * gsoft / (Eb + Eb0);
    double dG_T = -(Sb0 * delta / kappa) * (pow(1.0 + Eb / Eb0, bigN) * exp((Temp - Temp_0) / kappa));
    double pE_evt = -dG_Eb * (Sb / gsoft);
    double pE_Tvt = -dG_T * (Sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * chi * xi * Sb;
  
    /* Obtain the timestep */
    double dt = GetTimeStep();

    fXxii = theta * dt * hh * pE_ptau;
  }
  return fXxii;
}

double tevp2D::ComputeCtcon(void)
{
  /* compute Ctcon - implement the imperfection into the viscoplasticity
   * stress accumulate function and calculate the evolution function */
  //cout << "ComputeCtcon called" << endl;
  double Sb = fInternal[kSb];
  double Eb = fInternal[kEb];
  double Temp = fInternal[kTemp];
			      
  if (Sb <= kYieldTol)
    fCtcon = 0.0;
  else
  {
    dArrayT PP, Smlp;
    Smlp = GetSmlp();
    PP = GetPP();
    double pCp2 = dArrayT::Dot(Smlp, PP);
    pCp2 += pCp;
    double gsoft = Sb0 * pow(1.0 + Eb/Eb0, bigN);
    gsoft *= (1.0 - delta * (exp((Temp - Temp_0)/kappa) - 1.0));
    double reg = Sb / gsoft;
    double Ebtot = Eb0tot * pow(reg, smm);
    double pE_ptau = smm * Ebtot / Sb;
    double dG_Eb = bigN * gsoft / (Eb + Eb0);
    double dG_T = -(Sb0 * delta / kappa) * (pow(1.0 + Eb / Eb0, bigN) * exp((Temp - Temp_0) / kappa));
    double pE_evt = -dG_Eb * (Sb / gsoft);
    double pE_Tvt = -dG_T * (Sb / gsoft);
    double hh = pCp2 - pE_evt - pE_Tvt * chi * xi * Sb;
    
    /* Obtain the timestep */
    double dt = GetTimeStep();

    double xxii = theta * dt * hh * pE_ptau;
    fCtcon = xxii / ((1.0 + xxii) * hh);
  }
  return fCtcon;
}

dArrayT& tevp2D::ComputeSmlp(void)
{
  /* used in ComputePP */
  /* compute the deviatoric Kirchoff stress */
  //cout << "ComputeSmlp called" << endl;
  dMatrixT temp_stress(3); 
  temp_stress = ArrayToMatrix(fTempStress);
  CauchyToKirchoff(temp_stress);
  double trace_KH = fKirchoff.Trace() / 3.0;
  fSmlp[0] = fKirchoff(0,0) - trace_KH;
  fSmlp[1] = fKirchoff(1,1) - trace_KH;
  fSmlp[2] = fKirchoff(2,2) - trace_KH;
  fSmlp[3] = fKirchoff(0,1);   // Stored in Voigt notation
  double Sb = fInternal[kSb];
  if (Sb <= kYieldTol)
    int blah = 0;
  else
    fSmlp *= 1.5 / Sb;
  return fSmlp;
}

dArrayT& tevp2D::ComputePP(void)
{
  //cout << "ComputePP called" << endl;
  fPP = 0.0;
  dArrayT smlp(kVoigt);
  dMatrixT Dmat(kVoigt);
  smlp = GetSmlp();
  Dmat = GetDmat();
  Dmat.Multx(smlp, fPP);
  return fPP;
}

double tevp2D::ComputeEcc(void)
{
  /* Access ecc */
  //cout << "ComputeEcc called" << endl;
  dArrayT PP(kVoigt);
  PP = GetPP();
  ComputeGradients();
  if (fInternal[kSb] <= kYieldTol)   // If hasn't yielded yet...
    fEcc = 0.0;
  else
  { 
    for (int i = 0; i < 3; i++) 
      fEcc += PP[i] * fDtot(i,i);

    fEcc += PP[3] * fDtot(0,1);
  }
  return fEcc;
}

dMatrixT& tevp2D::ComputeDmat(void)
{
  /* computes the original elastic coefficient tensor */
  //cout << "ComputeDmat called" << endl;
  fDmat = 0.0;
  double Ed = el_E * (1.0 - el_V) / ((1.0 + el_V) * (1.0 - 2.0 * el_V));
  double Es = el_E * el_V / ((1.0 + el_V) * (1.0 - 2.0 * el_V));
  double G0 = el_E / (1.0 + el_V);
  fDmat(0,0) = fDmat(1,1) = fDmat(2,2) = Ed;
  fDmat(0,1) = fDmat(1,0) = fDmat(2,0) = fDmat(2,1) = Es;
  fDmat(0,2) = fDmat(0,3) = fDmat(1,2) = fDmat(1,3) = 0.0;
  fDmat(2,3) = fDmat(3,0) = fDmat(3,1) = fDmat(3,2) = 0.0;
  fDmat(3,3) = G0;

  return fDmat;
}

dArrayT& tevp2D::ComputeEP_tan(void)
{
  /* computes the modulus correction if plasticity has occurred */
  cout << "ComputeEP_tan called" << endl;
  dArrayT PP(kVoigt), diagU(kVoigt);
  double Ebtot = GetEbtot();
  double xxii = GetXxii();
  PP = GetPP();
  double Sb = fInternal[kSb];
  diagU[0] = 1.0;
  diagU[1] = 1.0;
  diagU[2] = 1.0;
  diagU[3] = 0.0;  

  for (int i = 0; i < 4; i++) {
    /* EP_tan is the plastic corrector to the tangent modulus */
    fEP_tan[i] = (Ebtot / (1.0 + xxii)) * (PP[i] + 3.0 * el_K * alpha_T * xi * chi * Sb * diagU[i]);  
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
  int TotalIP = fNumIP;                // Get total# of gauss points
  int i_size = 0;
  int d_size = 0;
  i_size += 2 * TotalIP;              // 2 flags per IP:  critical strain
                                      // and check for plasticity
  d_size += kNumOutput * TotalIP;     // 3 internal variables to track
  d_size += kVoigt * TotalIP;         // 4 non-zero stress components
                                      // Sig11, Sig12=Sig21, Sig22 and Sig33
  /* construct new plastic element */
  element.Allocate(i_size, d_size);

  /* first set of flags for plasticity criterion */
  for (int ip = 0; ip < TotalIP; ip++)
    (element.IntegerData())[ip] = kIsElastic;
  
  /* second set of flags for critical strain / model switch criterion */
  for (int ip = TotalIP; ip < 2 * TotalIP; ip++)
    (element.IntegerData())[ip] = kTevp;

  element.DoubleData() = 0.0;
}

void tevp2D::LoadData(const ElementCardT& element, int ip)
{
  /* load element data for the specified integration point */
  /* check */
  //cout << "LoadData called" << endl;
  int TotalIP = fNumIP;
  if (!element.IsAllocated()) throw eGeneralFail;

  int dex = ip * kVoigt;     // 4 non-zero stress components (11, 12, 22, 33)
  int offset = TotalIP * 4;

  /* fetch arrays */
  dArrayT& d_array = element.DoubleData();
  fTempStress.Set(kVoigt, &d_array[dex]);
  fInternal.Set(kNumOutput, &d_array[offset + ip * kNumOutput]);
}

void tevp2D::Update(ElementCardT& element)
{
  /* get flags */
  cout << "Update called" << endl;
  iArrayT& flags = element.IntegerData();
  int TotalIP = fNumIP;
  /* check if reset state (is same for all ip) */
  if (flags[0] == kReset)
  {
      flags = kIsElastic;          // don't update again
      return;
  }

  /* update plastic variables */
  for (int ip = 0; ip < TotalIP; ip++)
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

dArrayT tevp2D::MatrixToArray(dSymMatrixT StressMatrix)
{
  /* Flattens Kirchoff stress matrix into array form for internal variable
   * storage */
  //cout << "MatrixToArray called" << endl;
  fStressArray[0] = StressMatrix[0];
  fStressArray[3] = StressMatrix[5]; 
  fStressArray[1] = StressMatrix[1];
  fStressArray[2] = StressMatrix[2];
  
  return fStressArray;
}

dMatrixT tevp2D::ArrayToMatrix(dArrayT StressArray)
{
  /* Expands internal variable stress array to matrix form */
  //cout << "ArrayToMatrix called" << endl;
  fStressMatrix = 0.0;
  fStressMatrix(0,0) = StressArray[0];
  fStressMatrix(1,0) = fStressMatrix(0,1) = StressArray[3];
  fStressMatrix(1,1) = StressArray[1];
  fStressMatrix(2,2) = StressArray[2];
  return fStressMatrix;
}

dSymMatrixT tevp2D::Return3DStress(dMatrixT StressMatrix)
{
  /* Takes 3D matrix and converts to 3D symmetric matrix - necessary
   * because canned functions depend on fNumSD to convert */
  //cout << "Return3DStress called" << endl;
  fStill3D[0] = StressMatrix(0,0);
  fStill3D[1] = StressMatrix(1,1);
  fStill3D[2] = StressMatrix(2,2);
  fStill3D[3] = StressMatrix(1,2);
  fStill3D[4] = StressMatrix(0,2);
  fStill3D[5] = StressMatrix(0,1);
  return fStill3D;
}





