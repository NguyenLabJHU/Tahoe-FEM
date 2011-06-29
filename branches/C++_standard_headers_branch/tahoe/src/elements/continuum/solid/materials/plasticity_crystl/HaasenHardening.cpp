/*
  File: HaasenHardening.cpp
*/

#include "HaasenHardening.h"
#include "PolyCrystalMatT.h"
#include "NLCSolver.h"
#include "Utils.h"


// number of material properties and initial hard values

using namespace Tahoe;

const int kNumMatProp = 8;
const int kNumInitValues = 3;

HaasenHardening::HaasenHardening(PolyCrystalMatT& poly) :
  SlipHardening(poly)
{
  // input file
  ifstreamT& in = poly.Input_x();

  // number of hardening variables
  fNumHardVar = poly.NumSlip();

  // allocate space for ...
  // ... material properties and initial values of hard variables 
  fMatProp.Dimension(kNumMatProp);
  fInitHardValues.Dimension(kNumInitValues);

  // ... coeffs in hardening laws
  fAaa.Dimension(fNumHardVar);
  fC0aa.Dimension(fNumHardVar);
  fCaa.Dimension(fNumHardVar);
  fDa.Dimension(fNumHardVar);

  // ... hardening variables at t and t_n
  fDDtot.Dimension(fNumHardVar);
  fDDimm.Dimension(fNumHardVar);
  fXaa.Dimension(fNumHardVar);
  fDDtot_n.Dimension(fNumHardVar);
  fDDimm_n.Dimension(fNumHardVar);
  fXaa_n.Dimension(fNumHardVar);

  // ... accumulated slip shear strain at t and t_n
  fGamma.Dimension(fNumHardVar);
  fGamma_n.Dimension(fNumHardVar);

  // ... nondirectional hardness and effective stress
  fTauIso.Dimension(fNumHardVar);
  fTauEff.Dimension(fNumHardVar);

  // ... slip interaction and latent hardening matrices
  fFab.Dimension(fNumHardVar);
  fLatentMatx.Dimension(fNumHardVar);

  // ... workspace
  farray.Dimension(fNumHardVar);
 
  // input material properties for hardening law
  in >> fMatProp[0];   // Burger's vector "b" [m]
  in >> fMatProp[1];   // constant G (~ shear modulus) [Mpa]
  in >> fMatProp[2];   // constant K [m^(-1) MPa^(-1)]
  in >> fMatProp[3];   // latent hardening coeff "q" [adimensional]
  in >> fMatProp[4];   // gamma0 [adimensional]
  in >> fMatProp[5];   // X0aa [adimensional]
  in >> fMatProp[6];   // strain rate sensit exponent "m" [adimensional]
  in >> fMatProp[7];   // tau_0 [MPa]

  fMatProp[2] *= fMatProp[7];   // constant K [m^(-1)]

  // input initial values of slip system hardening variables
  in >> fInitHardValues;     // fDDtot_0, fDDimm_0, fXaa_0

  // compute temperature dependent coefficients
  TempDependentFuncs(in);

  // compute slip interaction matrix
  SlipInteractionMatrix(in);

  // build the latent hardening matrix
  LatentHardMatrix();

  // set hardening solver (in base class)
  SetHardeningSolver(in, 3*fNumHardVar);
}

HaasenHardening::~HaasenHardening() {}

void HaasenHardening::InitializeHardVariables() 
{ 
  // initialize hardening values at t_n
  fDDtot_n = fInitHardValues[0]; 
  fDDimm_n = fInitHardValues[1]; 
  fXaa_n   = fInitHardValues[2]; 

  // initialize total slip shear strain at t_n
  fGamma_n = 0.;

  // intialize variables at t
  ResetHistory();
}

double HaasenHardening::Magnitude() const { return fTauIso.Magnitude(); }
const int HaasenHardening::NumberOfVariables() const { return 9*fNumHardVar; }

void HaasenHardening::UpdateHistory() 
{ 
  // update hardening variables
  fDDtot_n = fDDtot; 
  fDDimm_n = fDDimm; 
  fXaa_n   = fXaa; 

  // update total slip shear strain
  fGamma_n = fGamma;
}

void HaasenHardening::ResetHistory() 
{
  // reset hardening variables
  fDDtot = fDDtot_n; 
  fDDimm = fDDimm_n; 
  fXaa   = fXaa_n; 

  // reset total slip shear strain
  fGamma = fGamma_n;

  // reset nondirectional hardening stress
  NonDirectionalHardening(fTauIso, fDDtot_n, fXaa_n);
}

void HaasenHardening::LoadHardData(int dim, int dex, dArrayT& d_array)
{
  // recover hardening variables for current element/IP/grain
  fDDtot_n.Set (fNumHardVar,  &d_array[dex += dim        ]);
  fDDimm_n.Set (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fXaa_n.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fDDtot.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fDDimm.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fXaa.Set     (fNumHardVar,  &d_array[dex += fNumHardVar]);

  // recover total slip shear strain (assumes: fNumHardVar = fNumSlip!!)
  fGamma_n.Set (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fGamma.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);

  // nondirectional hardening stress
  fTauIso.Set  (fNumHardVar,  &d_array[dex += fNumHardVar]);
}

void HaasenHardening::ExplicitUpdateHard()
{
  // preliminary computations
  InternalHardQntsConst();
  InternalHardQntsVaria(fDDtot_n, fXaa_n);

  // forward Euler estimate for hardening variables
  for (int i = 0; i < fNumHardVar; i++)
    {
      fDDtot[i] = fDDtot_n[i] + fdt * HardeningLawDDtot(fDDtot_n[i], i);
      fDDimm[i] = fDDimm_n[i] + fdt * HardeningLawDDimm(fDDimm_n[i], i);
      fXaa[i]   = fXaa_n[i]   + fdt * HardeningLawXaa(fXaa_n[i], i);
    }

  // initial norm of hardening arrays
  fNormDDt0 = fDDtot.Magnitude();
  fNormDDi0 = fDDimm.Magnitude();
  fNormXaa0 = fXaa.Magnitude();
}

void HaasenHardening::ImplicitUpdateHard()
{
  // preliminary computations
  InternalHardQntsConst();
  InternalHardQntsVaria(fDDtot, fXaa);

  // implicit update for hardening variables
  for (int i = 0; i < fNumHardVar; i++)
    {
      fDDtot[i] = fDDtot_n[i] + fdt * HardeningLawDDtot(fDDtot[i], i);
      fDDimm[i] = fDDimm_n[i] + fdt * HardeningLawDDimm(fDDimm[i], i);
      fXaa[i]   = fXaa_n[i]   + fdt * HardeningLawXaa(fXaa[i], i);
    }

  // norm of hardening arrays
  fNormDDt = fDDtot.Magnitude();
  fNormDDi = fDDimm.Magnitude();
  fNormXaa = fXaa.Magnitude();
}

void HaasenHardening::ImplicitSolveHard()
{
  // backward Euler estimate
  throwRunTimeError("HaasenHardening::ImplicitSolveHard called");
}

void HaasenHardening::FormRHS(const dArrayT& tauIso, dArrayT& rhs)
{
  // form residual
  #pragma unused(tauIso, rhs)
  throwRunTimeError("HaasenHardening::FormRHS called");
}

void HaasenHardening::FormLHS(const dArrayT& tauIso, dMatrixT& lhs)
{
  // form jacobian
  #pragma unused(tauIso, lhs)
  throwRunTimeError("HaasenHardening::FormLHS called");
}

bool HaasenHardening::Converged(double toler)
{
  // bool variable to check convergence on hardening variables
  bool test = ( fabs( fNormDDt/fNormDDt0-1.0 ) < toler &&
                fabs( fNormDDi/fNormDDi0-1.0 ) < toler &&
                fabs( fNormXaa/fNormXaa0-1.0 ) < toler ); 

  // reset norms if did not converge
  if (!test) {
     fNormDDt0 = fNormDDt;
     fNormDDi0 = fNormDDi;
     fNormXaa0 = fNormXaa;
  }

  return test;
} 

double HaasenHardening::HardeningModulus() const
{
  // for now return G
  return  fMatProp[1];
}

const dArrayT& HaasenHardening::ComputeHardQnts()
{
  // compute term in Haasen's kinetic equation
  farray.SetToCombination(1., fDDtot, -1., fDDimm);
  farray *= (fMatProp[0]*fBT);
  return farray;
}

const double HaasenHardening::ComputeHardQnts(int is)
{
  // compute term "is" in Haasen's kinetic equation
  return fMatProp[0]*fBT*(fDDtot[is]-fDDimm[is]);
}

/* PRIVATE MEMBERS FUNCTIONS */

void HaasenHardening::NonDirectionalHardening(dArrayT& tauIso, const dArrayT& DDtot,
					      const dArrayT& Xaa)
{
  // auxiliar quantity: array_i = Xab_i*Xab_i*DDtot_i
  for (int i = 0; i < fNumHardVar; i++)
    farray[i] = Xaa[i] * Xaa[i] * DDtot[i];

  // compute the non directional slip system hardness
  for (int i = 0; i < fNumHardVar; i++)
    tauIso[i] = fMatProp[0] * fMatProp[1] * sqrt(fLatentMatx.DotRow(i, farray));
}

void HaasenHardening::TempDependentFuncs(ifstreamT& in)
{
  // some constants
  double B0 = 9.50e-5, B1 = 2.67e3;  // [m s^-1 (MPa)^(-1/m)], [K]
  double B2 = 1.83e5,  B3 = 2.19e4;  // [m s^-1 (MPa)^(-1/m)], [K]
  double G0 = 5.45e5,  G1 = 2.19e4;  // adimensional, [K]
  double W0 = 1.09e7,  W1 = 2.19e4;  // [1/s], [K]

  // input temperature
  double theta;
  in >> theta;

  // compute functions B(T), G(T), and W(T)
  fBT = B0*exp(-B1/theta) + B2*exp(-B3/theta);
  fGT = G0*exp(-G1/theta);
  fWT = W0*exp(-W1/theta);

  fBT *= pow(fMatProp[7], 1./fMatProp[6]);  // [m s^(-1)]
}

void HaasenHardening::SlipInteractionMatrix(ifstreamT& in)
{
  // input corrector coeffs of interaction matrix
  in >> fC1xFab >> fC2xFab;

  // input Bassani-Wu coeffs of interaction matrix
  double r1, r2, r3, r4, r5;
  in >> r1 >> r4 >> r5;
  r2 = r3 = r1;

  // interaction matrix
  double* pF = fFab.Pointer();

  // ... column 0
  *pF++ = 0.; *pF++ = r3; *pF++ = r3; *pF++ = r5; *pF++ = r4; *pF++ = r2; 
     *pF++ = r5; *pF++ = r2; *pF++ = r4; *pF++ = r4; *pF++ = r4; *pF++ = r1;

  // ... column 1
  *pF++ = r3; *pF++ = 0.; *pF++ = r3; *pF++ = r2; *pF++ = r4; *pF++ = r5; 
     *pF++ = r4; *pF++ = r4; *pF++ = r1; *pF++ = r5; *pF++ = r2; *pF++ = r4;

  // ... column 2
  *pF++ = r3; *pF++ = r3; *pF++ = 0.; *pF++ = r4; *pF++ = r1; *pF++ = r4; 
     *pF++ = r2; *pF++ = r5; *pF++ = r4; *pF++ = r2; *pF++ = r5; *pF++ = r4;

  // ... column 3
  *pF++ = r5; *pF++ = r2; *pF++ = r4; *pF++ = 0.; *pF++ = r3; *pF++ = r3; 
     *pF++ = r5; *pF++ = r4; *pF++ = r2; *pF++ = r4; *pF++ = r1; *pF++ = r4;

  // ... column 4
  *pF++ = r4; *pF++ = r4; *pF++ = r1; *pF++ = r3; *pF++ = 0.; *pF++ = r3; 
     *pF++ = r2; *pF++ = r4; *pF++ = r5; *pF++ = r2; *pF++ = r4; *pF++ = r5;

  // ... column 5
  *pF++ = r2; *pF++ = r5; *pF++ = r4; *pF++ = r3; *pF++ = r3; *pF++ = 0.; 
     *pF++ = r4; *pF++ = r1; *pF++ = r4; *pF++ = r5; *pF++ = r4; *pF++ = r2;

  // ... column 6
  *pF++ = r5; *pF++ = r4; *pF++ = r2; *pF++ = r5; *pF++ = r2; *pF++ = r4; 
     *pF++ = 0.; *pF++ = r3; *pF++ = r3; *pF++ = r1; *pF++ = r4; *pF++ = r4;

  // ... column 7
  *pF++ = r2; *pF++ = r4; *pF++ = r5; *pF++ = r4; *pF++ = r4; *pF++ = r1; 
     *pF++ = r3; *pF++ = 0.; *pF++ = r3; *pF++ = r4; *pF++ = r5; *pF++ = r2;

  // ... column 8
  *pF++ = r4; *pF++ = r1; *pF++ = r4; *pF++ = r2; *pF++ = r5; *pF++ = r4; 
     *pF++ = r3; *pF++ = r3; *pF++ = 0.; *pF++ = r4; *pF++ = r2; *pF++ = r5;

  // ... column 9
  *pF++ = r4; *pF++ = r5; *pF++ = r2; *pF++ = r4; *pF++ = r2; *pF++ = r5; 
     *pF++ = r1; *pF++ = r4; *pF++ = r4; *pF++ = 0.; *pF++ = r3; *pF++ = r3;

  // ... column 10
  *pF++ = r4; *pF++ = r2; *pF++ = r5; *pF++ = r1; *pF++ = r4; *pF++ = r4; 
     *pF++ = r4; *pF++ = r5; *pF++ = r2; *pF++ = r3; *pF++ = 0.; *pF++ = r3;

  // ... column 11
  *pF++ = r1; *pF++ = r4; *pF++ = r4; *pF++ = r4; *pF++ = r5; *pF++ = r2; 
     *pF++ = r4; *pF++ = r2; *pF++ = r5; *pF++ = r3; *pF++ = r3; *pF = 0.;
}

void HaasenHardening::LatentHardMatrix()
{
  // build the latent hardening matrix
  fLatentMatx = (fMatProp[3]*fMatProp[3]);
  for (int i = 0; i < fNumHardVar; i++)
    fLatentMatx(i, i) = 1.;
}

void HaasenHardening::InternalHardQntsConst()
{
  // update total slip shear strain Gamma (Gamma >= 0)
  for (int i = 0; i < fNumHardVar; i++)
    fGamma[i] = fGamma_n[i] + fabs(fDGamma[i]);

  // coefficient fCaa/fC0aa in evolution eqn for Xaa
  for (int i = 0; i < fNumHardVar; i++)
    farray[i] = tanh(fGamma[i]/fMatProp[4]);

  for (int i = 0; i < fNumHardVar; i++)
    fCaa[i] = 1. + fC2xFab * fFab.DotRow(i, farray);
}

void HaasenHardening::InternalHardQntsVaria(dArrayT& DDtot, dArrayT& Xaa)
{
  // non-directional slip hardening
  NonDirectionalHardening(fTauIso, DDtot, Xaa);

  // effective stress and coefficient Aaa (Xaa's evolution eqn)
  for (int i = 0; i < fNumHardVar; i++)
    {
      fTauEff[i] = fabs(fTau[i]) - fTauIso[i];
      if (fTauEff[i] <= 0.) fTauEff[i] = 0.;
      fAaa[i] = 0.5 + pow(fTauEff[i]/5.0, 2.5);
    }

  // coefficient C0aa (Xaa's evolution eqn)
  for (int i = 0; i < fNumHardVar; i++)
    //fC0aa[i] = 20. + 0.05 * sinh(15.*(1.-Xaa[i]/fAaa[i]));
    fC0aa[i] = 20.;

  // coefficient Da (DDimm's evolution eqn)
  for (int i = 0; i < fNumHardVar; i++)
    fDa[i] = fC1xFab * fFab.DotRow(i, fDDtot);
}
	
const double HaasenHardening::HardeningLawDDtot(double DDtot, int is)
{
  return  fMatProp[2]*DDtot*fBT*pow(fTauEff[is]/fMatProp[7], 1./fMatProp[6]+1.);
}

const double HaasenHardening::HardeningLawDDimm(double DDimm, int is)
{
  return  fDa[is]*fabs(fDGamma[is])/fdt - fWT*DDimm;
}

const double HaasenHardening::HardeningLawXaa(double Xaa, int is)
{
  return  fC0aa[is]*fCaa[is]*(1.-Xaa/fAaa[is])*fabs(fDGamma[is])/fdt 
          - fGT*(Xaa-fMatProp[5]);
}
