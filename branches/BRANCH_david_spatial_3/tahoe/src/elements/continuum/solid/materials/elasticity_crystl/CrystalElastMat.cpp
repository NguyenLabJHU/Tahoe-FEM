/* $Id: CrystalElastMat.cpp,v 1.6 2004-07-15 08:28:16 paklein Exp $ */
/*
  File: CrystalElastMat.cpp
*/

#include "CrystalElastMat.h"
#include "CrystalElast.h"

#include <iostream.h>
#include <math.h>

#include "StringT.h"

#include "dMatrixT.h"
#include "ArrayT.h"
#include "dArrayT.h"

#include "Utils.h"


using namespace Tahoe;

CrystalElastMat::CrystalElastMat(CrystalElast& poly)
{
#pragma unused(poly)
  // initialize material coefficients
  fC11 = fC12 = fC44 = 0.;
}

CrystalElastMat::~CrystalElastMat() { }

void CrystalElastMat::ElasticityProps(dArrayT& matprop, double Temp_DegC, int elem, int intpt)
{
#pragma unused(elem)
#pragma unused(intpt)

  CalculateModuli(Temp_DegC);

  // material properties used in computations
  matprop[0] = fC44;                     // mu 
  matprop[1] = fC12;                     // lambda
  matprop[2] = 0.5*(2.*fC44+fC12-fC11);  // beta (=0 for isotropic)
}

void CrystalElastMat::ThermalProps(dMatrixT& alpha, double Temp_DegC)
{
  // get thermal expansion coefficients
  CalculateAlpha(alpha, Temp_DegC);
}

bool CrystalElastMat::IsIsotropic() const { return false; }

//Temperature dependent experimental expressions for solid CdZnTe
//from Queheillalt and Wadley (1998) with temperature range
//of 0 degrees Celsius and 1098 degrees Celsius (melting point)

void CrystalElastMat::CalculateModuli(double DegC)
{
  // check temperature range
  bool test = ( DegC <= 1098. && DegC >= 0.);
  if (!test) throwRunTimeError("CrystalElastMat::CalculateModuli: Bad values of DegC");
  // Temperature dependent elastic stiffness constants
  fC11 = TempDepModuli(DegC, 56.166, -8.497e-03, -5.292e-06);
  fC12 = TempDepModuli(DegC, 36.785, -6.715e-03, -2.229e-06);
  fC44 = TempDepModuli(DegC, 20.242, -3.182e-03, -5.429e-07);
}

void  CrystalElastMat::CalculateAlpha(dMatrixT& alpha, double DegC)
{
  double* palpha = alpha.Pointer();
  double DegK = DegC + 273.15;

  // check temperature range
  bool test = ( DegC <= 1098. && DegC >= 0.);
  if (!test) throwRunTimeError("CrystalElastMat::CalculateAlpha: Bad values of DegC");

  //Temperature dependent thermal expansion coefficients
  palpha[0] = 5.345e-06 + 8.373e-10*DegK;
  palpha[4] = palpha[8] = palpha[0];
}

double CrystalElastMat::TempDepModuli(double Temp, double const1, double const2
, double const3)
{
  double GPa = const1+const2*Temp+const3*Temp*Temp;
  return GPa*pow(10.0, 9);
}
