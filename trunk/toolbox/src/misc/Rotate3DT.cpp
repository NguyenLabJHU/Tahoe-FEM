/* This class provides the functionality to do 3D coordinate */
/* transformations.                                          */

/* ....TO BE CHECKED (SA 07/25/02)....*/

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
Rotate3DT::Rotate3DT(void): fAngleDeg1(0.0), fAngleDeg2(0.0),fAngle1(0.0), fAngle2(0.0),
			    fQ(kMatrixDim), fVec(kMatrixDim)
{

}

Rotate3DT::Rotate3DT(double angle1,double angle2): fAngleDeg1(angle1), fAngleDeg2(angle2),
						   fQ(kMatrixDim),fVec(kMatrixDim)
{
	SetAngle(fAngleDeg1,fAngleDeg2);
}

void Rotate3DT::SetAngle(double angle1,double angle2)
{
	fAngle1 = angle1*Pi/180.0;
	fAngle2 = angle2*Pi/180.0;

	fQ(0,0) = cos(angle1)*sin(angle2);
	fQ(1,0) = sin(angle1)*sin(angle2);
	fQ(2,0) = cos(angle2);

	fQ(0,1) = cos(angle1)*cos(angle2);
	fQ(1,1) = sin(angle1)*cos(angle2);
	fQ(2,1) = -sin(angle2);

	fQ(0,2) = fQ(1,0)*fQ(2,1)-fQ(2,0)*fQ(1,1);
	fQ(1,2) =-(fQ(0,0)*fQ(2,1)-fQ(2,0)*fQ(0,1));
	fQ(2,2) = fQ(0,0)*fQ(1,1)-fQ(1,0)*fQ(0,1);
}
