/*
  File: SlipGeometry.cpp
*/

#include "SlipGeometry.h"
#include "Utils.h"

using namespace Tahoe;

/* number of slip systems based on crystal structure */
const int kSlipFCC   = 12;
const int kSlipBCC   = 12;
const int kSlipHCP   = 12;
const int kSlipFCC24 = 24;   // differentiate (+) and (-) slip directions
const int kSlipPE    = 9;    //1=b, 2= c
const int kSlipPEa   = 9;   // 1=a, 2=c

/* Base Class */


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

FCCGeometry::FCCGeometry():
  SlipGeometry(kSlipFCC)
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
  if (fNumSlip != kSlipFCC)
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

BCCGeometry::BCCGeometry():
  SlipGeometry(kSlipBCC)
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
  if (fNumSlip != kSlipBCC)
    throwRunTimeError("BCCGeometry::SetSlipVectors: NumSlip != 12");

  // Miller indices normal to the slip plane 
  fVecM[0][0] = 0.; fVecM[0][1] = 1.; fVecM[0][2] = 1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 0.; fVecM[1][2] = 1.;
  fVecM[2][0] = 1.; fVecM[2][1] = -1.; fVecM[2][2] = 0.;

  fVecM[3][0] = 0.; fVecM[3][1] = 1.; fVecM[3][2] =-1.;
  fVecM[4][0] = 1.; fVecM[4][1] = 0.; fVecM[4][2] = 1.;
  fVecM[5][0] = 1.; fVecM[5][1] = 1.; fVecM[5][2] = 0.;

  fVecM[6][0] = 0.; fVecM[6][1] = 1.; fVecM[6][2] = 1.;
  fVecM[7][0] = 1.; fVecM[7][1] = 0.; fVecM[7][2] =-1.;
  fVecM[8][0] = 1.; fVecM[8][1] = 1.; fVecM[8][2] = 0.;

  fVecM[ 9][0] = 0.; fVecM[ 9][1] = 1.; fVecM[ 9][2] =-1.;
  fVecM[10][0] = 1.; fVecM[10][1] = 0.; fVecM[10][2] =-1.;
  fVecM[11][0] = 1.; fVecM[11][1] =-1.; fVecM[11][2] = 0.;

  // Miller indices tangent to the slip plane
  fVecS[0][0] = 1.; fVecS[0][1] = 1.; fVecS[0][2] =-1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 1.; fVecS[1][2] =-1.;
  fVecS[2][0] = 1.; fVecS[2][1] = 1.; fVecS[2][2] =-1.;

  fVecS[3][0] = 1.; fVecS[3][1] =-1.; fVecS[3][2] =-1.;
  fVecS[4][0] = 1.; fVecS[4][1] =-1.; fVecS[4][2] =-1.;
  fVecS[5][0] = 1.; fVecS[5][1] =-1.; fVecS[5][2] =-1.;

  fVecS[6][0] = 1.; fVecS[6][1] =-1.; fVecS[6][2] = 1.;
  fVecS[7][0] = 1.; fVecS[7][1] =-1.; fVecS[7][2] = 1.;
  fVecS[8][0] = 1.; fVecS[8][1] =-1.; fVecS[8][2] = 1.;

  fVecS[ 9][0] = 1.; fVecS[ 9][1] = 1.; fVecS[ 9][2] = 1.;
  fVecS[10][0] = 1.; fVecS[10][1] = 1.; fVecS[10][2] = 1.;
  fVecS[11][0] = 1.; fVecS[11][1] = 1.; fVecS[11][2] = 1.;
}


/* Derived Class HCPGeometry */

HCPGeometry::HCPGeometry():
  SlipGeometry(kSlipHCP)
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
  if (fNumSlip != kSlipHCP)
    throwRunTimeError("HCPGeometry::SetSlipVectors: NumSlip != 12");

  // ...
  fVecM[0][0] = 2.; fVecM[0][1] =-1.; fVecM[0][2] = 1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 2.; fVecM[1][2] =-1.;
  fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] = 2.;

  fVecM[3][0] = 2.; fVecM[3][1] = 1.; fVecM[3][2] = 1.;
  fVecM[4][0] = 1.; fVecM[4][1] = 2.; fVecM[4][2] =-1.;
  fVecM[5][0] = 1.; fVecM[5][1] =-1.; fVecM[5][2] = 2.;

  fVecM[6][0] = 2.; fVecM[6][1] = 1.; fVecM[6][2] =-1.;
  fVecM[7][0] = 1.; fVecM[7][1] = 2.; fVecM[7][2] = 1.;
  fVecM[8][0] = 1.; fVecM[8][1] =-1.; fVecM[8][2] =-2.;

  fVecM[ 9][0] = 2.; fVecM[ 9][1] =-1.; fVecM[ 9][2] =-1.;
  fVecM[10][0] = 1.; fVecM[10][1] =-2.; fVecM[10][2] = 1.;
  fVecM[11][0] = 1.; fVecM[11][1] = 1.; fVecM[11][2] =-2.;

  // ...
  fVecS[0][0] = 1.; fVecS[0][1] = 1.; fVecS[0][2] =-1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 1.; fVecS[1][2] =-1.;
  fVecS[2][0] = 1.; fVecS[2][1] = 1.; fVecS[2][2] =-1.;

  fVecS[3][0] = 1.; fVecS[3][1] =-1.; fVecS[3][2] =-1.;
  fVecS[4][0] = 1.; fVecS[4][1] =-1; fVecS[4][2]  =-1.;
  fVecS[5][0] = 1.; fVecS[5][1] =-1.; fVecS[5][2] =-1.;

  fVecS[6][0] = 1.; fVecS[6][1] =-1.; fVecS[6][2] = 1.;
  fVecS[7][0] = 1.; fVecS[7][1] =-1.; fVecS[7][2] = 1.;
  fVecS[8][0] = 1.; fVecS[8][1] =-1.; fVecS[8][2] = 1.;

  fVecS[ 9][0] = 1.; fVecS[ 9][1] = 1.; fVecS[ 9][2] = 1.;
  fVecS[10][0] = 1.; fVecS[10][1] = 1.; fVecS[10][2] = 1.;
  fVecS[11][0] = 1.; fVecS[11][1] = 1.; fVecS[11][2] = 1.;
}


/* Derived Class FCCGeometry24 */

// distinguishes between positive and negative slip directions
FCCGeometry24::FCCGeometry24():
  SlipGeometry(kSlipFCC24)
{
  // set Miller indices for FCC - 24 slip systems
  SetSlipVectors();

  // compute slip system quantities
  InitializeSlipQnts();

  /*  for (int i = 0; i < fNumSlip; i++)
    {
      cout << " fZ_cr, Slip System # " << i + 1 << endl;
      cout << fZ[i] << "\n\n";
      } */
}

FCCGeometry24::~FCCGeometry24() {}

void FCCGeometry24::PrintName(ostream& out) const
{
  //print crystal structure type
  out << "    FCC24 crystal structure\n";
}

void FCCGeometry24::SetSlipVectors()
{
  // check for bad input value
  if (fNumSlip != kSlipFCC24)
    throwRunTimeError("FCCGeometry24::SetSlipVectors: NumSlip != 24");

  // Miller indices normal to the slip plane
  fVecM[0][0] = 1.; fVecM[0][1] = 1.; fVecM[0][2] =-1.;
  fVecM[1][0] = 1.; fVecM[1][1] = 1.; fVecM[1][2] =-1.;
  fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] =-1.;
  fVecM[3][0] = 1.; fVecM[3][1] = 1.; fVecM[3][2] =-1.;
  fVecM[4][0] = 1.; fVecM[4][1] = 1.; fVecM[4][2] =-1.;
  fVecM[5][0] = 1.; fVecM[5][1] = 1.; fVecM[5][2] =-1.;

  fVecM[ 6][0] = 1.; fVecM[ 6][1] =-1.; fVecM[ 6][2] =-1.;
  fVecM[ 7][0] = 1.; fVecM[ 7][1] =-1.; fVecM[ 7][2] =-1.;
  fVecM[ 8][0] = 1.; fVecM[ 8][1] =-1.; fVecM[ 8][2] =-1.;
  fVecM[ 9][0] = 1.; fVecM[ 9][1] =-1.; fVecM[ 9][2] =-1.;
  fVecM[10][0] = 1.; fVecM[10][1] =-1.; fVecM[10][2] =-1.;
  fVecM[11][0] = 1.; fVecM[11][1] =-1.; fVecM[11][2] =-1.;

  fVecM[12][0] = 1.; fVecM[12][1] =-1.; fVecM[12][2] = 1.;
  fVecM[13][0] = 1.; fVecM[13][1] =-1.; fVecM[13][2] = 1.;
  fVecM[14][0] = 1.; fVecM[14][1] =-1.; fVecM[14][2] = 1.;
  fVecM[15][0] = 1.; fVecM[15][1] =-1.; fVecM[15][2] = 1.;
  fVecM[16][0] = 1.; fVecM[16][1] =-1.; fVecM[16][2] = 1.;
  fVecM[17][0] = 1.; fVecM[17][1] =-1.; fVecM[17][2] = 1.;

  fVecM[18][0] = 1.; fVecM[18][1] = 1.; fVecM[18][2] = 1.;
  fVecM[19][0] = 1.; fVecM[19][1] = 1.; fVecM[19][2] = 1.;
  fVecM[20][0] = 1.; fVecM[20][1] = 1.; fVecM[20][2] = 1.;
  fVecM[21][0] = 1.; fVecM[21][1] = 1.; fVecM[21][2] = 1.;
  fVecM[22][0] = 1.; fVecM[22][1] = 1.; fVecM[22][2] = 1.;
  fVecM[23][0] = 1.; fVecM[23][1] = 1.; fVecM[23][2] = 1.;

  // Miller indices tangent to the slip plane
  fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 1.;
  fVecS[1][0] = 1.; fVecS[1][1] = 0.; fVecS[1][2] = 1.;
  fVecS[2][0] = 1.; fVecS[2][1] =-1.; fVecS[2][2] = 0.;
  fVecS[3][0] = 0.; fVecS[3][1] =-1.; fVecS[3][2] =-1.;
  fVecS[4][0] =-1.; fVecS[4][1] = 0.; fVecS[4][2] =-1.;
  fVecS[5][0] =-1.; fVecS[5][1] = 1.; fVecS[5][2] = 0.;

  fVecS[ 6][0] = 0.; fVecS[ 6][1] = 1.; fVecS[ 6][2] =-1.;
  fVecS[ 7][0] = 1.; fVecS[ 7][1] = 0.; fVecS[ 7][2] = 1.;
  fVecS[ 8][0] = 1.; fVecS[ 8][1] = 1.; fVecS[ 8][2] = 0.;
  fVecS[ 9][0] = 0.; fVecS[ 9][1] =-1.; fVecS[ 9][2] = 1.;
  fVecS[10][0] =-1.; fVecS[10][1] = 0.; fVecS[10][2] =-1.;
  fVecS[11][0] =-1.; fVecS[11][1] =-1.; fVecS[11][2] = 0.;

  fVecS[12][0] = 0.; fVecS[12][1] = 1.; fVecS[12][2] = 1.;
  fVecS[13][0] = 1.; fVecS[13][1] = 0.; fVecS[13][2] =-1.;
  fVecS[14][0] = 1.; fVecS[14][1] = 1.; fVecS[14][2] = 0.;
  fVecS[15][0] = 0.; fVecS[15][1] =-1.; fVecS[15][2] =-1.;
  fVecS[16][0] =-1.; fVecS[16][1] = 0.; fVecS[16][2] = 1.;
  fVecS[17][0] =-1.; fVecS[17][1] =-1.; fVecS[17][2] = 0.;

  fVecS[18][0] = 0.; fVecS[18][1] = 1.; fVecS[18][2] =-1.;
  fVecS[19][0] = 1.; fVecS[19][1] = 0.; fVecS[19][2] =-1.;
  fVecS[20][0] = 1.; fVecS[20][1] =-1.; fVecS[20][2] = 0.;
  fVecS[21][0] = 0.; fVecS[21][1] =-1.; fVecS[21][2] = 1.;
  fVecS[22][0] =-1.; fVecS[22][1] = 0.; fVecS[22][2] = 1.;
  fVecS[23][0] =-1.; fVecS[23][1] = 1.; fVecS[23][2] = 0.;
}

/* Derived Class PEGeometry 
 according to Aurelie Azoug's notes*/
PEGeometry::PEGeometry():
SlipGeometry(kSlipPE)
{
    // set Miller indices for HCP
    SetSlipVectors();
    
    // compute slip system quantities
    InitializeSlipQnts();
}

PEGeometry::~PEGeometry() {}

void PEGeometry::PrintName(ostream& out) const
{
    //print crystal structure type
    out << "    Polyethylene crystal structure\n";
}

void PEGeometry::SetSlipVectors()
{
    // check for bad input value
    if (fNumSlip != kSlipPE)
        throwRunTimeError("PEGeometry::SetSlipVectors: NumSlip != 9");
    
    // ...permutated From Aurelie's slides
    fVecM[0][0] = 0.; fVecM[0][1] = 0.; fVecM[0][2] =1.;    /*A*/
    fVecM[1][0] = 1.; fVecM[1][1] = 0.; fVecM[1][2] =0.;    /*B*/ /*Crystal chain slip*/
    fVecM[2][0] = 1.; fVecM[2][1] = 0.; fVecM[2][2] =1.;    /*C*/
    
    fVecM[3][0] = 0.; fVecM[3][1] =0.; fVecM[3][2] =1.;     /*D*/
    fVecM[4][0] = 1.; fVecM[4][1] =0.; fVecM[4][2] =0.;     /*E*/   /*Crystal transverse slip*/
    fVecM[5][0] = 1.; fVecM[5][1] =0.; fVecM[5][2] =1.;     /*F*/
    
    fVecM[6][0] = 0.; fVecM[6][1] =1.; fVecM[6][2] = 0.;    /*K*/   /*Kink Bands*/
    fVecM[7][0] = 0.; fVecM[7][1] =1.; fVecM[7][2] = 0.;    /*K*/
    
    fVecM[8][0] = 0.; fVecM[8][1] =1.; fVecM[8][2] = 0.;    /*Chain Pull*/
    
    // ...
    fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 0.;
    fVecS[1][0] = 0.; fVecS[1][1] = 1.; fVecS[1][2] = 0.;
    fVecS[2][0] = 0.; fVecS[2][1] = 1.; fVecS[2][2] = 0.;
    
    fVecS[3][0] = 1.; fVecS[3][1] = 0.; fVecS[3][2] = 0.;
    fVecS[4][0] = 0.; fVecS[4][1] = 0.; fVecS[4][2] = 1.;
    fVecS[5][0] = -1.; fVecS[5][1] = 0.; fVecS[5][2] = 1.;
    
    fVecS[6][0] = 0.; fVecS[6][1] = 0.; fVecS[6][2] = 1.;
    fVecS[7][0] = 1.; fVecS[7][1] = 0.; fVecS[7][2] = 0.;
    
    fVecS[8][0] = 0.; fVecS[8][1] = 1.; fVecS[8][2] = 0.;

    
    
#if (0)
    // ...From Aurelie's slides
    fVecM[0][0] = 1.; fVecM[0][1] = 0.; fVecM[0][2] =0.;    /*A*/
    fVecM[1][0] = 0.; fVecM[1][1] = 1.; fVecM[1][2] =0.;    /*B*/ /*Crystal chain slip*/
    fVecM[2][0] = 1.; fVecM[2][1] = 1.; fVecM[2][2] =0.;    /*C*/
    
    fVecM[3][0] = 1.; fVecM[3][1] =0.; fVecM[3][2] =0.;     /*D*/
    fVecM[4][0] = 0.; fVecM[4][1] =1.; fVecM[4][2] =0.;     /*E*/   /*Crystal transverse slip*/
    fVecM[5][0] = 1.; fVecM[5][1] =1.; fVecM[5][2] =0.;     /*F*/
    
    fVecM[6][0] = 0.; fVecM[6][1] =0.; fVecM[6][2] = 1.;    /*K*/   /*Kink Bands*/
    fVecM[7][0] = 0.; fVecM[7][1] =0.; fVecM[7][2] = 1.;    /*K*/
    
    fVecM[8][0] = 0.; fVecM[8][1] =0.; fVecM[8][2] = 1.;    /*Chain Pull*/
    
    // ...
    fVecS[0][0] = 0.; fVecS[0][1] = 0.; fVecS[0][2] = 1.;
    fVecS[1][0] = 0.; fVecS[1][1] = 0.; fVecS[1][2] = 1.;
    fVecS[2][0] = 0.; fVecS[2][1] = 0.; fVecS[2][2] = 1.;
    
    fVecS[3][0] = 0.; fVecS[3][1] = 1.; fVecS[3][2] = 0.;
    fVecS[4][0] = 1.; fVecS[4][1] = 0.; fVecS[4][2] = 0.;
    fVecS[5][0] = 1.; fVecS[5][1] = -1.; fVecS[5][2] = 0.;
    
    fVecS[6][0] = 1.; fVecS[6][1] = 0.; fVecS[6][2] = 0.;
    fVecS[7][0] = 0.; fVecS[7][1] = 1.; fVecS[7][2] = 0.;
    
    fVecS[8][0] = 0.; fVecS[8][1] = 0.; fVecS[8][2] = 1.;
#endif
}

/* Derived Class PEGeometry
 according to Aurelie Azoug's notes*/
PEGeometryA::PEGeometryA():
SlipGeometry(kSlipPEa)
{
    // set Miller indices for HCP
    SetSlipVectors();
    
    // compute slip system quantities
    InitializeSlipQnts();
}

PEGeometryA::~PEGeometryA() {}

void PEGeometryA::PrintName(ostream& out) const
{
    //print crystal structure type
    out << "    Polyethylene crystal structure\n";
}

void PEGeometryA::SetSlipVectors()
{
    // check for bad input value
    if (fNumSlip != kSlipPEa)
        throwRunTimeError("PEGeometryA::SetSlipVectors: NumSlip != 9");
    
    // ...permutated From Aurelie's slides
    fVecM[0][0] = 1.; fVecM[0][1] = 0.; fVecM[0][2] = 0.; /*A*/
    fVecM[1][0] = 0.; fVecM[1][1] = 0.; fVecM[1][2] = 1.; /*B*/ /*Crystal chain slip*/
    fVecM[2][0] = 1.; fVecM[2][1] = 0.; fVecM[2][2] = 1.; /*C*/

    fVecM[3][0] = 1.; fVecM[3][1] = 0.; fVecM[3][2] = 0.; /*D*/
    fVecM[4][0] = 0.; fVecM[4][1] = 0.; fVecM[4][2] = 1.; /*E*/   /*Crystal transverse slip*/
    fVecM[5][0] = 1.; fVecM[5][1] = 0.; fVecM[5][2] = 1.; /*F*/
    
    fVecM[6][0] = 0.; fVecM[6][1] = 1.; fVecM[6][2] = 0.; /*K*/   /*Kink Bands*/
    fVecM[7][0] = 0.; fVecM[7][1] = 1.; fVecM[7][2] = 0.; /*K*/
    
    fVecM[8][0] = 0.; fVecM[8][1] = 1.; fVecM[8][2] = 0.; /*Chain Pull*/
    
    // ...
    fVecS[0][0] = 0.; fVecS[0][1] = 1.; fVecS[0][2] = 0.;
    fVecS[1][0] = 0.; fVecS[1][1] = 1.; fVecS[1][2] = 0.;
    fVecS[2][0] = 0.; fVecS[2][1] = 1.; fVecS[2][2] = 0.;
    
    fVecS[3][0] = 0.; fVecS[3][1] = 0.; fVecS[3][2] = 1.;
    fVecS[4][0] = 1.; fVecS[4][1] = 0.; fVecS[4][2] = 0.;
    fVecS[5][0] = 1.; fVecS[5][1] = 0.; fVecS[5][2] = -1.;
    fVecS[6][0] = 1.; fVecS[6][1] = 0.; fVecS[6][2] = 0.;
    
    fVecS[7][0] = 0.; fVecS[7][1] = 0.; fVecS[7][2] = 1.;
    fVecS[8][0] = 0.; fVecS[8][1] = 1.; fVecS[8][2] = 0.;

     
}


