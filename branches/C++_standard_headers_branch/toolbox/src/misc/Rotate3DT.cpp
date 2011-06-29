/* $Id: Rotate3DT.cpp,v 1.5 2003-11-10 22:14:34 cjkimme Exp $ */
/* This class provides the functionality to do 3D coordinate */
/* transformations.                                          */

#include "Rotate3DT.h"
#include <math.h>
#include "toolboxConstants.h"

#include "dArrayT.h"
#include "dArray2DT.h"

/* size parameters */

using namespace Tahoe;

const int kMatrixDim = 3;

const double Pi = acos(-1.0);

/*
* constructor
*/
Rotate3DT::Rotate3DT(void): fDir(kMatrixDim),fAngleDeg(0.0),fAngle(0.0), 
			    fQ(kMatrixDim), fVec(kMatrixDim),
			    fPhiDeg(0.0),fThetaDeg(0.0),fPsiDeg(0.0),
			    fPhi(0.0),fTheta(0.0),fPsi(0.0)
{

}

Rotate3DT::Rotate3DT(dArrayT u,double angle): fDir(kMatrixDim),fAngleDeg(angle), 
					      fAngle(0.0),fQ(kMatrixDim),fVec(kMatrixDim),
					      fPhiDeg(0.0), fThetaDeg(0.0),fPsiDeg(0.0),
					      fPhi(0.0),fTheta(0.0),fPsi(0.0)
{
  SetAngle(u,fAngleDeg);
}

Rotate3DT::Rotate3DT(double phi,double theta,double psi): fDir(kMatrixDim),fAngleDeg(0.0),
					    fAngle(0.0),fQ(kMatrixDim),fVec(kMatrixDim),
					    fPhiDeg(phi),fThetaDeg(theta),fPsiDeg(psi),					    
					    fPhi(0.0),fTheta(0.0),fPsi(0.0)
{
  SetEuler(fPhiDeg,fThetaDeg,fPsiDeg);
}

void Rotate3DT::SetAngle(dArrayT u,double angle)
{
  fAngle = angle*Pi/180.0; // angle in radians
  for (int i=0; i<kMatrixDim; i++) fDir[i] = u[i];

  fCos = cos(fAngle);
  fSin = sin(fAngle);

  
  for (int i=0; i<kMatrixDim; i++)
    {
      for (int j=0; j<kMatrixDim; j++)
	fQ(i,j) = (1.-fCos)*fDir[i]*fDir[j];
      fQ(i,i) += fCos;
    }
  
  
  fQ(1,0) += fSin*fDir[2];
  fQ(2,0) -= fSin*fDir[1];
  
  fQ(0,1) -= fSin*fDir[2];
  fQ(2,1) += fSin*fDir[0];
  
  fQ(0,2) += fSin*fDir[1];
  fQ(1,2) -= fSin*fDir[0];
}

void Rotate3DT::GiveTransfoMatrix(const dArray2DT Q)
{
  // A transformation matrix is given by the user.
  // Replace fQ in  SetAngle.
  for (int i=0; i<kMatrixDim; i++)
    for (int j=0; j<kMatrixDim; j++)
      fQ(i,j) = Q(i,j);  
}

void Rotate3DT::SetEuler(double phi,double theta,double psi)
{
  fPhi = phi*Pi/180.0; 
  fTheta = theta*Pi/180.0; 
  fPsi = psi*Pi/180.0; 

  fCosPhi= cos(fPhi);
  fSinPhi= sin(fPhi);
  fCosTheta= cos(fTheta);
  fSinTheta= sin(fTheta);
  fCosPsi= cos(fPsi);
  fSinPsi= sin(fPsi);


  fQ(0,0) = fCosPhi*fCosPsi - fSinPhi*fCosTheta*fSinPsi;
  fQ(1,0) =-fCosPhi*fSinPsi - fSinPhi*fCosTheta*fCosPsi;
  fQ(2,0) =-fSinPhi*fSinTheta;

  fQ(0,1) = fSinPhi*fCosPsi + fCosPhi*fCosTheta*fSinPsi;
  fQ(1,1) =-fSinPhi*fSinPsi + fCosPhi*fCosTheta*fCosPsi; 
  fQ(2,1) = fCosPhi*fSinTheta;

  fQ(0,2) = fSinTheta*fSinPsi;
  fQ(1,2) = fSinTheta*fCosPsi;
  fQ(2,2) = fCosTheta;
}
