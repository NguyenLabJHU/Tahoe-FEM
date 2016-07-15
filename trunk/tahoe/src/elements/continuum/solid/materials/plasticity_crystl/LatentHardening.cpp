/*
  File: LatentHardening.cpp
*/

#include "LatentHardening.h"
#include "PolyCrystalMatT.h"
#include "NLCSolver.h"
#include "Utils.h"


// number of material properties and initial hard values

using namespace Tahoe;

const int kNumMatProp    = 2; /*fMatProp[0] =b, fMatProp[1]=Q*/
int kNumInitValues = 1; /*tau^\alpha_0*/
int kNumInternal   = 2;

// codes for computing hardening law qnts
const int kFunc  = 0;  // HardFunc 
const int kdFunc = 1;  // d(HardFunc)/d(Crss)

// codes to access fInternal
const int kShearRate_n = 0;
const int kShearRate   = 1; /*\dot{p}*/

// some limits for total shear rate
const double SHR_MIN = 0.;
const double SHR_MAX = 1.0e6;

LatentHardening::LatentHardening(PolyCrystalMatT& poly) :
  SlipHardening(poly),
  fInternal (kNumInternal*poly.NumSlip()),
  farray    (poly.NumSlip())
{
  // input file
  ifstreamT& in = poly.Input_x();

  // number of hardening variables
   fNumHardVar = poly.NumSlip();
  //fNumHardVar = 1;

  // allocate space for
  // ... material constants and initial hard values
  fMatProp.Dimension(kNumMatProp);
  fInitHardValues.Dimension(kNumInitValues*fNumHardVar);

  // ... hardening variables
  fTauIso.Dimension(fNumHardVar);
  fTauIso_n.Dimension(fNumHardVar);
  ftauiso_save.Dimension(fNumHardVar);

  // input material properties for hardening law
  in >> fMatProp[0];   // b
  in >> fMatProp[1];   // Q

  // input initial slip system hardness (different for different system)
    for(int i=0; i<kNumInitValues*fNumHardVar; i++)
        in >> fInitHardValues[i];
/*
    cout << "\nnumval: "<<kNumInitValues*fNumHardVar;
    cout << "\nfInitHardValues: "<<fInitHardValues;
*/
  // set hardening solver (in base class)
  SetHardeningSolver(in, fNumHardVar);
}

LatentHardening::~LatentHardening() {}

void LatentHardening::InitializeHardVariables()
{
/*Setting \tau^\alpha and \tau^\alpha_n to inital strength values*/

// hardening stress
    for(int i=0; i<kNumInitValues*fNumHardVar; i++)
    {
        fTauIso_n[i] = fInitHardValues[i];
//       fTauIso   = fInitHardValues[0];
    }
    fTauIso = fTauIso_n;
// internal variables
    fInternal = 0.0;
}

double LatentHardening::Magnitude() const { return fTauIso.Magnitude(); }
const int LatentHardening::NumberOfVariables() const 
{
/*\tau^\alpha, \tau^\alpha_n, \dot{p}, \dot{p}_*, shear strength and effective flow rates*/
  return 2*fNumHardVar + kNumInternal; 
}

void LatentHardening::UpdateHistory() 
{ 
  fTauIso_n = fTauIso; 
  fInternal[kShearRate_n] = fInternal[kShearRate];
}

void LatentHardening::ResetHistory() { 
  fTauIso = fTauIso_n; 
  fInternal[kShearRate] = fInternal[kShearRate_n];
}

void LatentHardening::LoadHardData(int dim, int dex, dArrayT& d_array)
{
  // recover hardening variables for current element/IP/grain
  fTauIso_n.Set (fNumHardVar,  &d_array[dex += dim        ]);
  fTauIso.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fInternal.Set (kNumInternal*fNumHardVar, &d_array[dex += fNumHardVar]);
}

void LatentHardening::ExplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // forward Euler estimate for crss
  double c   = fdt * fMatProp[0];
  double g_n, g;
  for (int i = 0; i < fNumHardVar; i++)
  {
    /*get H_n*/
      g_n = fTauIso_n[i]/fInitHardValues[i];
      g = g_n + fdt*fMatProp[0]*(fMatProp[1] - g_n)*fInternal[kShearRate];
      fTauIso[i] = g * fInitHardValues[i];
  }
//    cout << "\nExplicit fTauIso: "<<fTauIso;

  // norm of explicit estimate
  fNormHard0 = fTauIso.Magnitude();
}

void LatentHardening::ImplicitUpdateHard()   // called from Algorithm 1
{
  // preliminary computations
  InternalHardQnts();
    
    /*backward Euler*/
    double g_n, g;
    for (int i = 0; i < fNumHardVar; i++)
    {
//        cout <<"\ngn: "<<g_n;
//        cout <<"\ninternal: "<<fInternal[kShearRate];
        /*get H_n*/
        g_n = fTauIso_n[i]/fInitHardValues[i];
        g = g_n + fdt*fMatProp[0]*fMatProp[1]*fInternal[kShearRate];
        g /= 1.0+fdt*fMatProp[0]*fInternal[kShearRate];
        fTauIso[i] = g * fInitHardValues[i];
    }
#if 0
    cout << "\nImplicit Solve Hard fTauIso: "<<fTauIso;
    cout << "\ng: "<<g;
    cout << "\nInternal shear rate: "<<fInternal[kShearRate];
#endif
    // norm of implicit estimate
    fNormHard = fTauIso.Magnitude();
//    cout << "\nImplicit fTauIso: "<<fTauIso;

/*OLD
    // generalized mid-point approximation for crss (theta=0.5)
  double c   = 0.5 * fdt * fMatProp[0];
//  double g_s = fTauIsoSat - fMatProp[1];
  double g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i];
     if ( (g_n/g_s) <= 1.0 )
     {
        g = g_n + c*( (1.0-g_n/g_s)*fInternal[kShearRate_n] + fInternal[kShearRate] );
        g /= ( 1.0 + c*fInternal[kShearRate]/g_s );
     }
     else
     {
        g = g_n;
     }

     fTauIso[i] = g + fMatProp[1];
  }
*/
}

void LatentHardening::ImplicitSolveHard()   // called from Algorithm 2
{
//  int ierr = 0;
  
  // preliminary computations
  InternalHardQnts();

  // backward Euler method: call solver method of NLCSolver class
  //fSolver->Solve(fSolverPtr, fTauIso, ierr);

  //if (ierr == 1)
  //  throwRunTimeError("LatentHardening::SolveImplicitHard: Convergence problems");

  // backward Euler approximation for crss
    /*backward Euler*/
    double g_n, g;
    for (int i = 0; i < fNumHardVar; i++)
    {
        /*get H_n*/
        g_n = fTauIso_n[i]/fInitHardValues[i];
        g = g_n + fdt*fMatProp[0]*fMatProp[1]*fInternal[kShearRate];
        g /= 1.0+fdt*fMatProp[0]*fInternal[kShearRate];
        fTauIso[i] = g * fInitHardValues[i];
    }

  // norm of implicit solution
  fNormHard = fTauIso.Magnitude();
}

void LatentHardening::FormRHS(const dArrayT& tauIso, dArrayT& rhs)
{
  // form residual
  for (int i = 0; i < fNumHardVar; i++)
    rhs[i] = tauIso[i] - fTauIso_n[i] - fdt * fInitHardValues[i]*HardeningLaw(tauIso[i]/fInitHardValues[i], kFunc);
}

void LatentHardening::FormLHS(const dArrayT& tauIso, dMatrixT& lhs)
{
  // form jacobian
  lhs = 0.0;
  for (int i = 0; i < fNumHardVar; i++)
    lhs(i,i) = 1. - fdt * fInitHardValues[i]*HardeningLaw(tauIso[i]/fInitHardValues[i], kdFunc);
}

bool LatentHardening::Converged(double toler)
{
  // check convergence on hardening variables
//    cout << "\nfNormHard: "<< fNormHard;
//    cout << "\nfNormHard0: "<< fNormHard0;
 //   cout << "\ndiff norm: "<<fNormHard-fNormHard0;
  bool test = ( fabs(fNormHard-fNormHard0) < toler*fInitHardValues[0] );

  // if did not converge, reset norm
  if (!test) fNormHard0 = fNormHard;

  return test;
} 

void LatentHardening::SaveCurrentSolution() { ftauiso_save = fTauIso; }
void LatentHardening::RestoreSavedSolution() { fTauIso = ftauiso_save; }

double LatentHardening::HardeningModulus() const
{
  // for now return h0
  return  fMatProp[0];
}

const double LatentHardening::IsoHardeningStress(int is) const
{
  // one isotropic hardening variable per grain
  #pragma unused(is)
  return fTauIso[is];
}

#if 0 //removing this for now. The parent class returns tauIso
const dArrayT& LatentHardening::ComputeHardQnts()
{
  // preliminary computations
  InternalHardQnts();

    //compute dHardLaw/dTauIso
    double dHdIso = fdt * HardeningLaw(fTauIso[0], kdFunc);

  // compute factor of dHardLaw/dDGamma
//  double dHdDGam = fInitHardValues[i]*HardeningLaw(fTauIso[0]/fInitHardValues[i], kFunc) / (fInternal[kShearRate]*fdt);

  // compute array (1-dHardLaw/dTauIso)^(-1)*dHardLaw/dDGamma
  for (int i = 0; i < fDGamma.Length(); i++)
    {
        double dHdDGam = fInitHardValues[i]*HardeningLaw(fTauIso[i]/fInitHardValues[i], kFunc) / (fInternal[kShearRate]*fdt);
      farray[i] = fdt * dHdDGam * fabs(fDGamma[i]) / fDGamma[i];
//      farray[i] *= 1./ (1. - dHdIso);
    }
  return farray;
}
#endif

#if 0
  // print hardening parameters
  out << "       Hardening rate (h0) . . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       Initial hardness (g0) . . . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       Saturation hardness (gs0) . . . . . . . . = " << fMatProp[2] << "\n";
  out << "       Saturation shear rate (gams0) . . . . . . = " << fMatProp[3] << "\n";
  out << "       Saturation exponent (xms) . . . . . . . . = " << fMatProp[4] << "\n";
#endif

/* PRIVATE MEMBERS FUNCTIONS */

void LatentHardening::InternalHardQnts()
{
  // accumulated shear rate
  fInternal[kShearRate] = 0.;
  for (int i = 0; i < fDGamma.Length(); i++) fInternal[kShearRate] += fabs(fDGamma[i]);
  fInternal[kShearRate] /= fdt; 

  // check limits on fShearRate
  if (fInternal[kShearRate] <= SHR_MIN) fInternal[kShearRate] = SHR_MIN;
  if (fInternal[kShearRate] >= SHR_MAX) fInternal[kShearRate] = SHR_MAX;

  // hardening saturation level
//  fTauIsoSat = fMatProp[2] * pow(fInternal[kShearRate]/fMatProp[3], fMatProp[4]);
}
	
const double LatentHardening::HardeningLaw(double H, int kcode)
{
//  double hard = fMatProp[0] / (fTauIsoSat - fMatProp[1]) * fInternal[kShearRate];
    double hard = fMatProp[0]*fInternal[kShearRate];
  switch(kcode)
    {
    case kFunc:
//      hard *= (fTauIsoSat - tauIso);
        hard *= (fMatProp[1]-H);
      break;
    case kdFunc:
      hard *= -1.0;
      break;
    default:
      throwRunTimeError("LatentHardening::HardeningLaw: Bad kcode");
    }

  return hard;
}
