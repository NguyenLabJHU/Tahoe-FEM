/*
  File: SimplePowerLaw.cpp
*/

#include "SimplePowerLaw.h"
#include "EVPFDBaseT.h"


using namespace Tahoe;

const int kNumMatProp = 5;

SimplePowerLaw::SimplePowerLaw(EVPFDBaseT& model)
{
  // fetch input file
  ifstreamT& in = model.Input();

  // allocate space for matl's property array
  fMatProp.Dimension(kNumMatProp);

  // read material properties for kinetic equation
  in >> fMatProp[0];   // s0
  in >> fMatProp[1];   // e0
  in >> fMatProp[2];   // edot0
  in >> fMatProp[3];   // n (hardening exp)
  in >> fMatProp[4];   // m (strain rate sensitivity exp)
}

SimplePowerLaw::~SimplePowerLaw() { }

double SimplePowerLaw::g(double eqp)
{
  return  fMatProp[0]*exp(fMatProp[3]*log(1.+eqp/fMatProp[1]));
}

double SimplePowerLaw::f(double sigma, double eqp)
{
  //  double g = fMatProp[0] * pow((eqp/fMatProp[1]+1.), fMatProp[3]);
  //  return  fMatProp[2] * (pow((sigma/g), 1./fMatProp[4]) - 1.);

  double syield0 = g(eqp);
  return fMatProp[2] * (exp(1./fMatProp[4]*log(sigma/syield0)) - 1.);
}

double SimplePowerLaw::DfDsigma(double sigma, double eqp)
{
  double eqpdot = f(sigma, eqp);
  return  (eqpdot+fMatProp[2])/(fMatProp[4]*sigma);
}

double SimplePowerLaw::DfDs(double sigma, double eqp)
{
  double eqpdot = f(sigma, eqp);
  return  -fMatProp[3]/fMatProp[4]*(eqpdot+fMatProp[2])/(eqp+fMatProp[1]);
}

double SimplePowerLaw::h(double eqpdot, double kappa)
{
#pragma unused(eqpdot)
#pragma unused(kappa)

  throw ExceptionT::kGeneralFail;
}

double SimplePowerLaw::DhDeqpdot(double eqpdot, double kappa)
{
#pragma unused(eqpdot)
#pragma unused(kappa)

  throw ExceptionT::kGeneralFail;
}

double SimplePowerLaw::DhDs(double eqpdot, double kappa)
{
#pragma unused(eqpdot)
#pragma unused(kappa)

  throw ExceptionT::kGeneralFail;
}

void SimplePowerLaw::Print(ostream& out) const
{
  // print kinetics equation data
  out << "       reference stress (s0) . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       reference strain (e0) . . . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       reference strain rate (edot0) . . . . . . = " << fMatProp[2] << "\n";
  out << "       hardening exponent (n). . . . . . . . . . = " << fMatProp[3] << "\n";
  out << "       rate sensitivity exponent (m) . . . . . . = " << fMatProp[4] << "\n";
}

void SimplePowerLaw::PrintName(ostream& out) const
{
  // output kinetic equation model name
  out << "    Simple Power Law Kinetic Equation\n";
}  
   
