/*
  File: SlipGeometry.cpp
*/

#include "SlipGeometry.h"
#include "Utils.h"

/* Base Class */


using namespace Tahoe;

SlipGeometry::SlipGeometry(int numslip) :
  fNumSlip(numslip)
{
  // allocate space for arrays
  fVecS.Dimension(fNumSlip);
  fVecM.Dimension(fNumSlip);
  fZ.Dimension(fNumSlip);
  for (int i = 0; i < fNumSlip; i++)
    {
      fVecS[i].Dimension(3);
      fVecM[i].Dimension(3);
      fZ[i].Dimension(3,3);
    }
}

SlipGeometry::~SlipGeometry() {}

void SlipGeometry::InitializeSlipQnts()
{
  // build slip systems unit vectors
  NormalizeSlipVectors();

  // build Schmidt tensor fVecS(x)fVecM
  SetSchmidtTensor();
}

void SlipGeometry::NormalizeSlipVectors()
{
  for(int i = 0; i < fNumSlip; i++)
    {
     fVecS[i].UnitVector();
     fVecM[i].UnitVector();
    }    
}

void SlipGeometry::SetSchmidtTensor()
{
  for(int i = 0; i < fNumSlip; i++)
    {
      fZ[i].Outer(fVecS[i], fVecM[i]);
    }
}


/* Derived Class FCCGeometry */

FCCGeometry::FCCGeometry(int numslip):
  SlipGeometry(numslip) 
{
  // set Miller indices for FCC
  SetSlipVectors();

  // compute slip system quantities 
  InitializeSlipQnts();

  /*  for (int i = 0; i < fNumSlip; i++)
    {
      cout << " fZ_cr, Slip System # " << i + 1 << endl;
      cout << fZ[i] << "\n\n";
      } */
}

FCCGeometry::~FCCGeometry() {}

void FCCGeometry::PrintName(ostream& out) const
{
  //print crystal structure type
  out << "    FCC crystal structure\n";
}

void FCCGeometry::SetSlipVectors()
{
  // check for bad input value
  if (fNumSlip != 12) 
    throwRunTimeError("FCCGeometry::SetSlipVectors: NumSlip != 12");

  // Miller indices normal to the slip plane 
  fVecM[0][0] = 1.; fVecM[0][1] = 1.; fVecM[0][2] =-1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 1.; fVecM[1][2] =-1.;
  fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] =-1.;

  fVecM[3][0] = 1.; fVecM[3][1] =-1.; fVecM[3][2] =-1.;
  fVecM[4][0] = 1.; fVecM[4][1] =-1.; fVecM[4][2] =-1.;
  fVecM[5][0] = 1.; fVecM[5][1] =-1.; fVecM[5][2] =-1.;

  fVecM[6][0] = 1.; fVecM[6][1] =-1.; fVecM[6][2] = 1.;
  fVecM[7][0] = 1.; fVecM[7][1] =-1.; fVecM[7][2] = 1.;
  fVecM[8][0] = 1.; fVecM[8][1] =-1.; fVecM[8][2] = 1.;

  fVecM[ 9][0] = 1.; fVecM[ 9][1] = 1.; fVecM[ 9][2] = 1.;
  fVecM[10][0] = 1.; fVecM[10][1] = 1.; fVecM[10][2] = 1.;
  fVecM[11][0] = 1.; fVecM[11][1] = 1.; fVecM[11][2] = 1.;

  // Miller indices tangent to the slip plane 
  fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 0.; fVecS[1][2] = 1.;
  fVecS[2][0] = 1.; fVecS[2][1] =-1.; fVecS[2][2] = 0.;

  fVecS[3][0] = 0.; fVecS[3][1] = 1.; fVecS[3][2] =-1.;
  fVecS[4][0] = 1.; fVecS[4][1] = 0.; fVecS[4][2] = 1.;
  fVecS[5][0] = 1.; fVecS[5][1] = 1.; fVecS[5][2] = 0.;

  fVecS[6][0] = 0.; fVecS[6][1] = 1.; fVecS[6][2] = 1.;
  fVecS[7][0] = 1.; fVecS[7][1] = 0.; fVecS[7][2] =-1.;
  fVecS[8][0] = 1.; fVecS[8][1] = 1.; fVecS[8][2] = 0.;

  fVecS[ 9][0] = 0.; fVecS[ 9][1] = 1.; fVecS[ 9][2] =-1.;
  fVecS[10][0] = 1.; fVecS[10][1] = 0.; fVecS[10][2] =-1.;
  fVecS[11][0] = 1.; fVecS[11][1] =-1.; fVecS[11][2] = 0.;
}


/* Derived Class BCCGeometry */

BCCGeometry::BCCGeometry(int numslip):
  SlipGeometry(numslip)
{
  // set Miller indices for BCC
  SetSlipVectors();

  // compute slip system quantities 
  InitializeSlipQnts();
}

BCCGeometry::~BCCGeometry() {}

void BCCGeometry::PrintName(ostream& out) const
{
  //print crystal structure type
  out << "    BCC crystal structure\n";
}

void BCCGeometry::SetSlipVectors()
{
  // check for bad input value
  if (fNumSlip != 12) 
    throwRunTimeError("BCCGeometry::SetSlipVectors: NumSlip != 12");

  // Miller indices normal to the slip plane 
  fVecM[0][0] = 1.; fVecM[0][1] = 1.; fVecM[0][2] =-1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 1.; fVecM[1][2] =-1.;
  fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] =-1.;

  fVecM[3][0] = 1.; fVecM[3][1] =-1.; fVecM[3][2] =-1.;
  fVecM[4][0] = 1.; fVecM[4][1] =-1.; fVecM[4][2] =-1.;
  fVecM[5][0] = 1.; fVecM[5][1] =-1.; fVecM[5][2] =-1.;

  fVecM[6][0] = 1.; fVecM[6][1] =-1.; fVecM[6][2] = 1.;
  fVecM[7][0] = 1.; fVecM[7][1] =-1.; fVecM[7][2] = 1.;
  fVecM[8][0] = 1.; fVecM[8][1] =-1.; fVecM[8][2] = 1.;

  fVecM[ 9][0] = 1.; fVecM[ 9][1] = 1.; fVecM[ 9][2] = 1.;
  fVecM[10][0] = 1.; fVecM[10][1] = 1.; fVecM[10][2] = 1.;
  fVecM[11][0] = 1.; fVecM[11][1] = 1.; fVecM[11][2] = 1.;

  // Miller indices tangent to the slip plane
  fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 0.; fVecS[1][2] = 1.;
  fVecS[2][0] = 1.; fVecS[2][1] =-1.; fVecS[2][2] = 0.;

  fVecS[3][0] = 0.; fVecS[3][1] = 1.; fVecS[3][2] =-1.;
  fVecS[4][0] = 1.; fVecS[4][1] = 0.; fVecS[4][2] = 1.;
  fVecS[5][0] = 1.; fVecS[5][1] = 1.; fVecS[5][2] = 0.;

  fVecS[6][0] = 0.; fVecS[6][1] = 1.; fVecS[6][2] = 1.;
  fVecS[7][0] = 1.; fVecS[7][1] = 0.; fVecS[7][2] =-1.;
  fVecS[8][0] = 1.; fVecS[8][1] = 1.; fVecS[8][2] = 0.;

  fVecS[ 9][0] = 1.; fVecS[ 9][1] = 1.; fVecS[ 9][2] = 1.;
  fVecS[10][0] = 1.; fVecS[10][1] = 1.; fVecS[10][2] = 1.;
  fVecS[11][0] = 1.; fVecS[11][1] = 1.; fVecS[11][2] = 1.;
}


/* Derived Class HCPGeometry */

HCPGeometry::HCPGeometry(int numslip):
  SlipGeometry(numslip)
{
  // set Miller indices for HCP
  SetSlipVectors();

  // compute slip system quantities 
  InitializeSlipQnts();
}

HCPGeometry::~HCPGeometry() {}

void HCPGeometry::PrintName(ostream& out) const
{
  //print crystal structure type
  out << "    HCP crystal structure\n";
}

void HCPGeometry::SetSlipVectors()
{
  // check for bad input value
  if (fNumSlip != 12)
    throwRunTimeError("HCPGeometry::SetSlipVectors: NumSlip != 12");

  // ...
  fVecM[0][0] = 1.; fVecM[0][1] = 1.; fVecM[0][2] =-1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 1.; fVecM[1][2] =-1.;
  fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] =-1.;

  fVecM[3][0] = 1.; fVecM[3][1] =-1.; fVecM[3][2] =-1.;
  fVecM[4][0] = 1.; fVecM[4][1] =-1.; fVecM[4][2] =-1.;
  fVecM[5][0] = 1.; fVecM[5][1] =-1.; fVecM[5][2] =-1.;

  fVecM[6][0] = 1.; fVecM[6][1] =-1.; fVecM[6][2] = 1.;
  fVecM[7][0] = 1.; fVecM[7][1] =-1.; fVecM[7][2] = 1.;
  fVecM[8][0] = 1.; fVecM[8][1] =-1.; fVecM[8][2] = 1.;

  fVecM[ 9][0] = 1.; fVecM[ 9][1] = 1.; fVecM[ 9][2] = 1.;
  fVecM[10][0] = 1.; fVecM[10][1] = 1.; fVecM[10][2] = 1.;
  fVecM[11][0] = 1.; fVecM[11][1] = 1.; fVecM[11][2] = 1.;

  // ...
  fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 0.; fVecS[1][2] = 1.;
  fVecS[2][0] = 1.; fVecS[2][1] =-1.; fVecS[2][2] = 0.;

  fVecS[3][0] = 0.; fVecS[3][1] = 1.; fVecS[3][2] =-1.;
  fVecS[4][0] = 1.; fVecS[4][1] = 0.; fVecS[4][2] = 1.;
  fVecS[5][0] = 1.; fVecS[5][1] = 1.; fVecS[5][2] = 0.;

  fVecS[6][0] = 0.; fVecS[6][1] = 1.; fVecS[6][2] = 1.;
  fVecS[7][0] = 1.; fVecS[7][1] = 0.; fVecS[7][2] =-1.;
  fVecS[8][0] = 1.; fVecS[8][1] = 1.; fVecS[8][2] = 0.;

  fVecS[ 9][0] = 1.; fVecS[ 9][1] = 1.; fVecS[ 9][2] = 1.;
  fVecS[10][0] = 1.; fVecS[10][1] = 1.; fVecS[10][2] = 1.;
  fVecS[11][0] = 1.; fVecS[11][1] = 1.; fVecS[11][2] = 1.;
}
