/*
  File: SlipKinetics.cpp
*/

#include "SlipKinetics.h"
#include "PolyCrystalMatT.h"

SlipKinetics::SlipKinetics(PolyCrystalMatT& poly):
  fHard (poly.GetSlipHardening())
{

}

SlipKinetics::~SlipKinetics()
{

}
