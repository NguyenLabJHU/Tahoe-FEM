/*
  File: BCJKineticEqn.cpp
*/

#include "BCJKineticEqn.h"
#include "EVPFDBaseT.h"
#include "ifstreamT.h"


using namespace Tahoe;

const int kNumMatProp = 3;

BCJKineticEqn::BCJKineticEqn(EVPFDBaseT& model)
{
  // fetch input file
  ifstreamT& in = model.Input();

  // Temperature
  fTheta = model.Temperature();

  // allocate space for matl's property array
  fMatProp.Dimension(kNumMatProp);

  // read material constants for kinetic equation
  in >> fC1 >> fC2 >> fC3 >> fC4 >> fC5 >> fC6;
  in >> fC19 >> fC20 >> fC21;

  // material parameters
  ComputeMaterialProperties(fTheta); 
}

BCJKineticEqn::~BCJKineticEqn() { }

double BCJKineticEqn::g(double kappa)
{
  // static yield stress
  return  kappa+fMatProp[1];
}

double BCJKineticEqn::f(double sigma, double kappa)
{
  // kinetic equation function f
  double arg = sigma - kappa - fMatProp[1];
  if (arg > 0.0) 
     return  fMatProp[2]*sinh((sigma-kappa-fMatProp[1])/fMatProp[0]);
  else
     return 0.0;
}

double BCJKineticEqn::DfDsigma(double sigma, double kappa)
{
  // d(f)/d(s_eff)
  double arg = sigma - kappa - fMatProp[1];
  if (arg > 0.0) 
     return  fMatProp[2]/fMatProp[0]*cosh((sigma-kappa-fMatProp[1])/fMatProp[0]);
  else
     return 0.0;
}

double BCJKineticEqn::DfDs(double sigma, double kappa)
{
  // d(f)/d(kappa)
  double arg = sigma - kappa - fMatProp[1];
  if (arg > 0.0) 
     return  -fMatProp[2]/fMatProp[0]*cosh((sigma-kappa-fMatProp[1])/fMatProp[0]);
  else
     return 0.0;
}

double BCJKineticEqn::h(double eqpdot, double kappa)
{
  // dynamic yield stress h
  double arg = eqpdot/fMatProp[2];
  return  kappa + fMatProp[1] + fMatProp[0]*log(arg+sqrt(arg*arg+1.));
}

double BCJKineticEqn::DhDeqpdot(double eqpdot, double kappa)
{
#pragma unused(kappa)

  // d(h)/d(eqpdot)
  double arg = eqpdot/fMatProp[2];
  return  fMatProp[0]/fMatProp[2]/sqrt(arg*arg+1.) ;
}

double BCJKineticEqn::DhDs(double eqpdot, double kappa)
{
#pragma unused(eqpdot)
#pragma unused(kappa)

  // d(h)/d(kappa)
  return  1.0;
}

void BCJKineticEqn::Print(ostream& out) const
{
  // print temperature and material constants
  out << "       Temperature [K] . . . . . . . . . . . . . = " << fTheta << "\n";
  out << "       C1  [MPa] . . . . . . . . . . . . . . . . = " << fC1    << "\n";
  out << "       C2  [K] . . . . . . . . . . . . . . . . . = " << fC2    << "\n";
  out << "       C3  [MPa] . . . . . . . . . . . . . . . . = " << fC3    << "\n";
  out << "       C4  [K] . . . . . . . . . . . . . . . . . = " << fC4    << "\n";
  out << "       C5  [1/s] . . . . . . . . . . . . . . . . = " << fC5    << "\n";
  out << "       C6  [K] . . . . . . . . . . . . . . . . . = " << fC6    << "\n";
  out << "       C19 [1/K] . . . . . . . . . . . . . . . . = " << fC19   << "\n";
  out << "       C20 [K] . . . . . . . . . . . . . . . . . = " << fC20   << "\n";
  out << "       C21 [0] . . . . . . . . . . . . . . . . . = " << fC21   << "\n";

  // print material properties
  out << "       V   [MPa] . . . . . . . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       Y   [MPa] . . . . . . . . . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       f   [1/s] . . . . . . . . . . . . . . . . = " << fMatProp[2] << "\n";
}

void BCJKineticEqn::PrintName(ostream& out) const
{
  // output kinetic equation model name
  out << "    BCJ Kinetic Equation\n";
}
  
void BCJKineticEqn::ComputeMaterialProperties(const double theta)
{
  // compute material parameters for kinetic equation
  fMatProp[0]  = fC1*exp(-fC2/theta);               // V
  fMatProp[1]  = fC3/(fC21+exp(-fC4/theta));        // Y
  if (fC19 > 0.0) 
      fMatProp[1] *= 0.5*(1.+tanh(fC19*(fC20-theta)));  // Y
  fMatProp[2]  = fC5*exp(-fC6/theta);               // f
}   
