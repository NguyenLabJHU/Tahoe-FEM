/* This class provides the functionality to do 3D coordinate */
/* transformations.                                          */

#include "Rotate3DT.h"
#include <math.h>
#include "Constants.h"

/* size parameters */

using namespace Tahoe;

const int kMatrixDim = 3;

const double Pi = acos(-1.0);

/*
* constructor
*/
Rotate3DT::Rotate3DT(void): fAngleDeg(0.0), fDir(kMatrixDim),fAngle(0.0), 
			    fQ(kMatrixDim), fVec(kMatrixDim)
{

}

Rotate3DT::Rotate3DT(dArrayT u,double angle): fAngleDeg(angle), fDir(kMatrixDim), 
					      fAngle(0.0),fQ(kMatrixDim),fVec(kMatrixDim)
{
  SetAngle(u,fAngleDeg);
}

void Rotate3DT::SetAngle(dArrayT u,double angle)
{
  fAngle = angle*Pi/180.0; // angle in radians
  for (int i=0; i<kMatrixDim; i++) fDir[i] = u[i];


  
  for (int i=0; i<kMatrixDim; i++)
    {
      for (int j=0; j<kMatrixDim; j++)
	fQ(i,j) = (1.-fAngle)*fDir[i]*fDir[j];
      fQ(i,i) += cos(fAngle);
    }
  
  
  fQ(1,0) += sin(fAngle)*fDir[2];
  fQ(2,0) -= sin(fAngle)*fDir[1];
  
  fQ(0,1) -= sin(fAngle)*fDir[2];
  fQ(2,1) += sin(fAngle)*fDir[0];
  
  fQ(0,2) += sin(fAngle)*fDir[1];
  fQ(1,2) -= sin(fAngle)*fDir[0];
}
