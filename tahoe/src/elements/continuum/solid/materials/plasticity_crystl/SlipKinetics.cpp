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

void SlipKinetics::SetUpRateSensitivity()
{  }

void SlipKinetics::ComputeRateSensitivity()
{  }

bool SlipKinetics::IsMaxRateSensitivity()
{  return true; }

void SlipKinetics::RestoreRateSensitivity()
{  }

