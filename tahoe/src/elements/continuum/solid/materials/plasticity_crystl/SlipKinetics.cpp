/*
  File: SlipKinetics.cpp
*/

#include "SlipKinetics.h"
#include "PolyCrystalMatT.h"


using namespace Tahoe;

SlipKinetics::SlipKinetics(PolyCrystalMatT& poly):
  fHard (poly.GetSlipHardening())
{

}

SlipKinetics::~SlipKinetics()
{

}

void SlipKinetics::MaxMinArgPowerLaw(const double &xm)
{
   double EXPON = 300.0;

   fArgMax = exp( EXPON/(1./xm-1.0)*log(10.0));
   fArgMin = exp(-EXPON/(1./xm-1.0)*log(10.0));
}

/*  evaluates  x^y */
double SlipKinetics::Power(const double &x, const double &y)
{
   double power;

   if (x == 0.0) {
      if (y > 0.0) 
         power = 0.e0;
      else if (y < 0.0)
         power = 1.e+300;
      else 
         power = 1.e0; 
   } 
   else {
      power = y * log10(fabs(x));
      if (power > 300.0)
         power = 1.e+300;
      else
         power = pow(10.e0, power);
      if (x < 0.0) power *= -1.0; 
   }

   return power;
} 

void SlipKinetics::CheckArgumentRange(double &arg, const double &sign)
{
   if (fabs(arg) < fArgMin) 
      arg = fArgMin*sign;
   else if (fabs(arg) > fArgMax) 
      arg = fArgMax*sign;
}

void SlipKinetics::SetUpRateSensitivity()
{  }

void SlipKinetics::ComputeRateSensitivity()
{  }

bool SlipKinetics::IsMaxRateSensitivity()
{  return true; }

void SlipKinetics::RestoreRateSensitivity()
{  }

