/* $Id: Rotate2DT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (07/21/1996)                                          */
/* This class provides the functionality to do 2D coordinate              */
/* transformations.                                                       */

#ifndef _ROTATE_T_H_
#define _ROTATE_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dArrayT.h"

class Rotate2DT
{
public:

	/* constructor - angle(degrees) */
	Rotate2DT(void);
	Rotate2DT(double angle);

	/* set the angle field and associated work variables - angle
	 * passed in degrees */
	void SetAngle(double angle);
	double Angle(void) const; /* angle returned in degrees */
	
	/* transformation tensor */
	const dMatrixT& Q(void) const;
	
	/* Transformations */
	
	/* vectors */
	const dArrayT& RotateVectorIn(const dArrayT& vector);
	const dArrayT& RotateVectorOut(const dArrayT& vector);
	
	/* reduced index symmetric matrices */
	const dSymMatrixT& RotateRedMatIn(const dSymMatrixT& redmat);
	const dSymMatrixT& RotateRedMatOut(const dSymMatrixT& redmat);

	/* matrices */
	const dMatrixT& RotateMatrixIn(const dMatrixT& matrix);
	const dMatrixT& RotateMatrixOut(const dMatrixT& matrix);

	/* push an index in */
	void PushIndexIn(dMatrixT& matrix, int index);
	
	/* reduced index 4th order tensors */
	void RotateRedTensorIn(dMatrixT& matrix);
	void RotateRedTensorOut(dMatrixT& matrix);
		
private:

	double		fAngleDeg;/* angle in degrees */
	double		fAngle;	  /* angle in radians */

	double		fCost;	/* cos(1*fAngle) */
	double      fSint;	/* sin(1*fAngle) */
	double		fCos2t;	/* cos(2*fAngle) */
	double      fSin2t; /* sin(2*fAngle) */
	double		fCos4t;	/* cos(4*fAngle) */
	double      fSin4t; /* sin(4*fAngle) */

	dMatrixT	fQ; 	/* rotation tensor */
	
	/* return values */
	dArrayT		fVec;
	dSymMatrixT	fRedMat;	
	dMatrixT	fMat;
	
private:

	/*
	 * Crunch the numbers - rotate direction set the direction either
	 * (+) or (-) fAngle for (1) and (-1), respectively.
	 */
	void TransformO42D(dMatrixT& matrix, int RotateDirection = 1);
};

/* inline functions */

/* transformation tensor */
inline const dMatrixT& Rotate2DT::Q(void) const { return (fQ); }

/* angle returned in degrees */
inline double Rotate2DT::Angle(void) const { return(fAngleDeg); }

/* vectors */
inline const dArrayT& Rotate2DT::RotateVectorIn(const dArrayT& vector)
{
	fQ.MultTx(vector, fVec);
	return(fVec);
}

inline const dArrayT& Rotate2DT::RotateVectorOut(const dArrayT& vector)
{
	fQ.Multx(vector, fVec);
	return(fVec);
}

/* matrices */
inline const dMatrixT& Rotate2DT::RotateMatrixIn(const dMatrixT& matrix)
{
	fMat.MultQTBQ(fQ,matrix);	
	return(fMat);
}

inline const dMatrixT& Rotate2DT::RotateMatrixOut(const dMatrixT& matrix)
{
	fMat.MultQBQT(fQ,matrix);	
	return(fMat);
}

#endif /* _ROTATE_T_H_ */
