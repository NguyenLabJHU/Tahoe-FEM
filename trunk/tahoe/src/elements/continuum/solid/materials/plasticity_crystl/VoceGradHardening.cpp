/*
  File: VoceGradHardening.cpp
*/

#include "VoceGradHardening.h"
#include "PolyCrystalMatT.h"
#include "NLCSolver.h"
#include "Utils.h"


// number of material properties for hard model 

using namespace Tahoe;

const int kNumMatProp    = 9;
const int kNumInitValues = 1;
const int kNumInternal   = 4;

// codes for computing hardening law qnts
const int kFunc  = 0;  // HardFunc 
const int kdFunc = 1;  // d(HardFunc)/d(tauS)

// codes to access array fInternal
const int kShearRate_n = 0;
const int kShearRate   = 1;
const int kWorkRate_n  = 2;
const int kWorkRate    = 3;

/* to debug */
const bool XTAL_MESSAGES = false;

VoceGradHardening::VoceGradHardening(PolyCrystalMatT& poly):
  SlipHardening(poly),
  fInternal (kNumInternal)
{
  // input file
  ifstreamT& in = poly.Input_x();

  // number of hardening variables
  //fNumHardVar = poly.NumSlip();
  fNumHardVar = 1;

  // allocate space for
  // ... material constants and initial hard values
  fMatProp.Dimension(kNumMatProp);
  fInitHardValues.Dimension(kNumInitValues);

  // ... hardening variables
  fTauIso.Dimension(fNumHardVar);
  //fTauKin.Dimension(fNumHardVar);
  fTauIso_n.Dimension(fNumHardVar);
  fTauKin.Dimension(poly.NumSlip());
  fTauInc.Dimension(poly.NumSlip());

  // input material properties for hardening law
  in >> fMatProp[0];   // h0 
  in >> fMatProp[1];   // tauS0
  in >> fMatProp[2];   // tauSsat
  in >> fMatProp[3];   // c_s
  in >> fMatProp[4];   // c_x
  in >> fMatProp[5];   // lentgh scale l
  in >> fMatProp[6];   // burger's vector b
  in >> fMatProp[7];   // c_2
  // in >> fMatProp[8];   // c_1 (computed here)

  // reference to elastic material properties
  const dArrayT& elasProps = poly.MaterialProperties();

  // material constant c_1
  fMatProp[8] = 2.*fMatProp[0]/(fMatProp[3]*elasProps[0]*fMatProp[6]);

  // input initial value of statistically stored dislocation
  in >> fInitHardValues[0];

  // set hardening solver (in base class)
  SetHardeningSolver(in, fNumHardVar);
}

VoceGradHardening::~VoceGradHardening() {}

void VoceGradHardening::InitializeHardVariables() 
{ 
  // initialize hardening variables at t_0
  fTauIso_n = fInitHardValues[0];
  fTauIso   = fInitHardValues[0];

  // initialize internals
  fInternal = 0.0;
}

double VoceGradHardening::Magnitude() const { return fTauIso.Magnitude(); }
const int VoceGradHardening::NumberOfVariables() const 
{ 
  return 2*fNumHardVar + kNumInternal; 
}

void VoceGradHardening::UpdateHistory()
{
  // update hardening variables
  fTauIso_n = fTauIso;
  //fTauKin_n = fTauKin;  // keep track of fXe (in GradCrystalPlast)

  // update internals
  fInternal[kShearRate_n] = fInternal[kShearRate];
  fInternal[kWorkRate_n]  = fInternal[kWorkRate];
}

void VoceGradHardening::ResetHistory()
{
  // reset hardening variables
  fTauIso = fTauIso_n;
  //fTauKin = fTauKin_n;  // keep track of fXe (in GradCrystalPlast)

  // reset internals
  fInternal[kShearRate] = fInternal[kShearRate_n];
  fInternal[kWorkRate]  = fInternal[kWorkRate_n];
}

void VoceGradHardening::LoadHardData(int dim, int dex, dArrayT& d_array)
{
  // recover hardening/internal variables for current element/IP/grain
  fTauIso_n.Set (fNumHardVar,  &d_array[dex += dim        ]);
  fTauIso.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fInternal.Set (kNumInternal, &d_array[dex += fNumHardVar]);
}

void VoceGradHardening::ExplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // forward Euler estimate for fTauS
  //double scale = 50.e3*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //scale /= fMatProp[0];
  //double scale = 100.0*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //double beta = 3.125e4;
  //double beta = fMatProp[7]/fMatProp[8]*fMatProp[6]/fMatProp[5];
  //double scale = beta*fMatProp[3]/fMatProp[4];
  double scale = fMatProp[3]*fMatProp[7]/fMatProp[8];

  double c   = fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1] + 1.e-8;
     g = g_n + c*(1.0-g_n/g_s)*fInternal[kShearRate] 
             + c*scale*fInternal[kWorkRate]/g_n; 
     fTauIso[i] = g + fMatProp[1];
  }
}

void VoceGradHardening::ImplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // mid-point rule estimate for fTauIso
  //double scale = 50.e3*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //scale /= fMatProp[0];
  //double scale = 100.0*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //double beta = 3.125e4;
  //double beta = fMatProp[7]/fMatProp[8]*fMatProp[6]/fMatProp[5];
  //double scale = beta*fMatProp[3]/fMatProp[4];
  double scale = fMatProp[3]*fMatProp[7]/fMatProp[8];

  double c   = 0.5 * fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double A, B, C, g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1] + 1.e-8;
     
     A = 1.0 + c * fInternal[kShearRate] / g_s;
     B = - g_n - c * ( (1.0-g_n/g_s)*fInternal[kShearRate_n] 
                       + scale*fInternal[kWorkRate_n]/g_n 
                       + fInternal[kShearRate] );
     C = -c*scale*fInternal[kWorkRate]; 
     g = (- B + sqrt(B*B - 4.0*A*C)) / (2.*A);

     fTauIso[i] = g + fMatProp[1];
  }
}

void VoceGradHardening::ImplicitSolveHard()
{
  int ierr = 0;
  
  // preliminary computations
  InternalHardQnts();

  // compute DDss
  //fSolver->Solve(fSolverPtr, fTauIso, ierr);

  //if (ierr == 1)
  //  throwRunTimeError("VoceGradHardening::SolveImplicitHard: Convergence problems");

  // backward Euler estimate for fTauIso
  //double scale = 50.e3*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //scale /= fMatProp[0];
  //double scale = 100.0*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //double beta = 3.125e4;
  //double beta = fMatProp[7]/fMatProp[8]*fMatProp[6]/fMatProp[5];
  //double scale = beta*fMatProp[3]/fMatProp[4];
  double scale = fMatProp[3]*fMatProp[7]/fMatProp[8];

  double c   = fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double A, B, C, g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1] + 1.e-8;
     
     A = 1.0 + c * fInternal[kShearRate] / g_s;
     B = - g_n - c * fInternal[kShearRate];
     C = - c * scale * fInternal[kWorkRate]; 
     g = (- B + sqrt(B*B - 4.0*A*C)) / (2.*A);

     fTauIso[i] = g + fMatProp[1];
  }
}

void VoceGradHardening::FormRHS(const dArrayT& tauIso, dArrayT& rhs)
{
  // form residual
  for (int i = 0; i < fNumHardVar; i++)
    rhs[i] = tauIso[i] - fTauIso_n[i] - fdt * HardeningLaw(tauIso[i], kFunc);
}

void VoceGradHardening::FormLHS(const dArrayT& tauIso, dMatrixT& lhs)
{
  // form Jacobian
  lhs = 0.0;
  for (int i = 0; i < fNumHardVar; i++)
    lhs(i,i) = 1. - fdt * HardeningLaw(tauIso[i], kdFunc);
}

bool VoceGradHardening::Converged(double toler)
{
  #pragma unused(toler)
  return false;
}

double VoceGradHardening::HardeningModulus() const
{
  // temporary value
  return  fMatProp[0];
}

const double VoceGradHardening::IsoHardeningStress(int is) const
{
  // one isotropic hardening variable per grain
  #pragma unused(is)
  return fTauIso[0];
}

#if 0
  // print hardening parameters
  out << "       Hardening rate (h0) . . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       Initial hardness (tauS0). . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       Saturation hardness (tauSsat) . . . . . . = " << fMatProp[2] << "\n";
  out << "       c_s . . . . . . . . . . . . . . . . . . . = " << fMatProp[3] << "\n";
  out << "       c_x . . . . . . . . . . . . . . . . . . . = " << fMatProp[4] << "\n";
  out << "       length scale (l). . . . . . . . . . . . . = " << fMatProp[5] << "\n";
  out << "       burger's vector (b) . . . . . . . . . . . = " << fMatProp[6] << "\n";
  out << "       c_2 . . . . . . . . . . . . . . . . . . . = " << fMatProp[7] << "\n";
  out << "       c_1 (computed). . . . . . . . . . . . . . = " << fMatProp[8] << "\n";
#endif

/* PRIVATE MEMBER FUNCTIONS */

void VoceGradHardening::InternalHardQnts()
{
  // accumulated shear rate and work rate
  fInternal[kShearRate] = 0.;
  fInternal[kWorkRate]  = 0.;
  for (int i = 0; i < fDGamma.Length(); i++)
    {
      fInternal[kShearRate] += fabs(fDGamma[i]);
      fInternal[kWorkRate]  += fabs(fTauInc[i])*fabs(fDGamma[i]);
    }
  fInternal[kShearRate] /= fdt;
  fInternal[kWorkRate]  /= fdt;

  // saturation value of hardness
  fTauIsoSat = fMatProp[2];
}

const double VoceGradHardening::HardeningLaw(double tauIso, int kcode)
{
  double coeff = 50.e3*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //double coeff = 100.0*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  //double coeff = fMatProp[0]*1.0*fMatProp[3]/fMatProp[4]; // h0*Beta*c_s/c_x

  double diff  = tauIso - fMatProp[1] + 1.e-8; 
  double hard = fMatProp[0] / (fTauIsoSat - fMatProp[1]) * fInternal[kShearRate];

  switch(kcode)
    {
    case kFunc:
      hard *= (fTauIsoSat - tauIso);
      if (diff > 0.) hard += coeff/diff*fInternal[kWorkRate];
      break;

    case kdFunc:
      hard *= -1.0;
      if (diff > 0.) hard -= coeff/(diff*diff)*fInternal[kWorkRate];
      break;

    default:
      throwRunTimeError("VoceGradHardening::HardeningLaw: Bad kcode");
    }

  return hard;
}
