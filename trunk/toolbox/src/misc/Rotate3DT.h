/* This class provides the functionality to do 3D coordinate */
/* transformations. Similar to Rotate2D but less routines*/

/* Rotation in 3D defined by a direction  and an angle */

#ifndef _ROTATE3D_T_H_
#define _ROTATE3D_T_H_

/* direct members */
#include "dMatrixT.h"


#include "dArrayT.h"
#include "dArray2DT.h"


namespace Tahoe {

class Rotate3DT
{
public:

	/* constructor - angle(degrees) */
	Rotate3DT(void);
	Rotate3DT(dArrayT u,double angle);
	Rotate3DT(double phi, double theta, double psi);

	/* set the angle field and associated work variables - angle
	 * passed in degrees */
	void SetAngle(dArrayT u,double angle);
	double Angle(void) const; /* angle returned in degrees */
	dArrayT Direction(void) const;
	
	/* transformation tensor */
	const dMatrixT& Q(void) const;
	void GiveTransfoMatrix(const dArray2DT Q);
	void SetEuler(double phi,double theta,double psi);
	
	/* Transformations */
	
	/* vectors */
	const dArrayT& RotateVectorIn(const dArrayT& vector);
	const dArrayT& RotateVectorOut(const dArrayT& vector);
	
private:

	dArrayT fDir;
	double	fAngleDeg; /* angle in degrees */
	double	fAngle; /* angles in radians */

	dMatrixT	fQ; 
	dArrayT		fVec;

	double fPhiDeg,fThetaDeg,fPsiDeg;
	double fPhi,fTheta,fPsi;

	double fCos;
	double fSin;

	double fCosPhi,fSinPhi,fCosTheta,fSinTheta,fCosPsi,fSinPsi;
	
};

/* inline functions */

/* transformation tensor */
inline const dMatrixT& Rotate3DT::Q(void) const { return (fQ); }

/* angle returned in degrees */
inline double Rotate3DT::Angle(void) const { return(fAngleDeg); }
inline dArrayT Rotate3DT::Direction(void) const { return(fDir); }

/* vectors */
inline const dArrayT& Rotate3DT::RotateVectorIn(const dArrayT& vector)
{
	fQ.MultTx(vector, fVec);
	return(fVec);
}

inline const dArrayT& Rotate3DT::RotateVectorOut(const dArrayT& vector)
{
	fQ.Multx(vector, fVec);
	return(fVec);
}


} // namespace Tahoe 
#endif /* _ROTATE3D_T_H_ */
