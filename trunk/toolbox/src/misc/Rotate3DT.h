/* This class provides the functionality to do 3D coordinate */
/* transformations. Similar to Rotate2D but less routines*/



/* ....TO BE CHECKED (SA 07/25/02)....*/
 
#ifndef _ROTATE3D_T_H_
#define _ROTATE3D_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dArrayT.h"


namespace Tahoe {

class Rotate3DT
{
public:

	/* constructor - angle(degrees) */
	Rotate3DT(void);
	Rotate3DT(double angle1,double angle2);

	/* set the angle field and associated work variables - angle
	 * passed in degrees */
	void SetAngle(double angle1,double angle2);
	double Angle1(void) const; /* angle1 returned in degrees */
	double Angle2(void) const; /* angle2 returned in degrees */
	
	/* transformation tensor */
	const dMatrixT& Q(void) const;
	
	/* Transformations */
	
	/* vectors */
	const dArrayT& RotateVectorIn(const dArrayT& vector);
	const dArrayT& RotateVectorOut(const dArrayT& vector);
	
private:

	double	fAngleDeg1,fAngleDeg2; /* angles in degrees */
	double	fAngle1,fAngle2; /* angles in radians */

	dMatrixT	fQ; 
	dArrayT		fVec;
	
};

/* inline functions */

/* transformation tensor */
inline const dMatrixT& Rotate3DT::Q(void) const { return (fQ); }

/* angle returned in degrees */
inline double Rotate3DT::Angle1(void) const { return(fAngleDeg1); }
inline double Rotate3DT::Angle2(void) const { return(fAngleDeg2); }

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
