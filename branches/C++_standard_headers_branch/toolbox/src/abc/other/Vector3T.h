/* $Id: Vector3T.h,v 1.11 2005-07-29 08:11:58 paklein Exp $ */
/* created: paklein (02/11/2000) */
#ifndef _VECTOR_3_T_H_
#define _VECTOR_3_T_H_

/* Environmental */
#include "Environment.h"

namespace Tahoe {

/** utility class for 3D vector functions.
 * \note some functions do create temporary nTYPE instances */
template <class nTYPE>
class Vector3T
{
  public:

	/** default constructor */
	Vector3T(void);

	/** "copy" constructor. Copy the first 3 values from source */
	Vector3T(const nTYPE* source);

	/** \name accessors */
	/*@{*/
	/** element accessor */
	nTYPE& operator[](int dex);

	/** const element accessor */
	const nTYPE& operator[](int dex) const;

	/** type conversion operator. Convert a Vector3T to a (const nTYPE*) */
	operator const nTYPE*() const;

	/** type conversion operator. Convert a Vector3T to a (nTYPE*) */
	operator nTYPE*();
	/*@}*/

	/** \name assignment operator */
	/*@{*/
	Vector3T<nTYPE>& operator=(const nTYPE* rhs);
	Vector3T<nTYPE>& operator=(const nTYPE& value);
	/*@}*/

	/** \name mathematical operators */
	/*@{*/
	Vector3T<nTYPE>& operator+=(const nTYPE* rhs);
	Vector3T<nTYPE>& operator-=(const nTYPE* rhs);
	/*@}*/

	/** \name with scalars */
	/*@{*/
	Vector3T<nTYPE>& operator+=(const nTYPE& value);
	Vector3T<nTYPE>& operator-=(const nTYPE& value);
	Vector3T<nTYPE>& operator*=(const nTYPE& value);
	Vector3T<nTYPE>& operator/=(const nTYPE& value);
	/*@}*/

	/** \name average */
	/*@{*/
	Vector3T<nTYPE>& Average(const nTYPE* a, const nTYPE* b);
	Vector3T<nTYPE>& Average(const nTYPE* a, const nTYPE* b, const nTYPE* c);
	/*@}*/
	
	/** difference */
	Vector3T<nTYPE>& Diff(const nTYPE* a, const nTYPE* b);

	/** vector cross product */
	Vector3T<nTYPE>& Cross(const nTYPE* a, const nTYPE* b);

	/** inner product */
	static nTYPE Dot(const nTYPE* a, const nTYPE* b);

	/** L2 norm (between a and b) */
	nTYPE Norm(void) const;
	static nTYPE Norm(const nTYPE* a, const nTYPE* b);

	/** linear combinations <- a*A + b*B */
	Vector3T<nTYPE>& Combine(const nTYPE& a, const nTYPE* A,
	                         const nTYPE& b, const nTYPE* B);

	/** fill vector with random numbers in the range [-1 1]
	 * \param seed random number seed */
	void Random(int seed);

private:

	/** statically allocated data */
	nTYPE v[3];
};

/* inlines */
template <class nTYPE>
inline Vector3T<nTYPE>::Vector3T(void)
{ 
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
}

template <class nTYPE>
inline Vector3T<nTYPE>::Vector3T(const nTYPE* source)
{
	operator=(source);
}

/* accessor */
template <class nTYPE>
inline nTYPE& Vector3T<nTYPE>::operator[](int dex)
{
#if __option(extended_errorcheck)
	if (dex < 0 || dex > 2) ExceptionT::OutOfRange("Vector3T");
#endif
	return v[dex];
}

template <class nTYPE>
inline const nTYPE& Vector3T<nTYPE>::operator[](int dex) const
{
#if __option(extended_errorcheck)
	if (dex < 0 || dex > 2) ExceptionT::OutOfRange("Vector3T");
#endif
	return v[dex];
}

/* type conversion operators */
template <class nTYPE>
inline Vector3T<nTYPE>::operator const nTYPE*() const { return (const nTYPE*) v; }

template <class nTYPE>
inline Vector3T<nTYPE>::operator nTYPE*() { return (nTYPE*) v; }

/* assignment operator */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator=(const nTYPE* rhs)
{
	v[0] = rhs[0];
	v[1] = rhs[1];
	v[2] = rhs[2];
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator=(const nTYPE& value)
{
	v[0] = value;
	v[1] = value;
	v[2] = value;
	return *this;
}

/* mathematical operators */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator+=(const nTYPE* rhs)
{
	v[0] += rhs[0];
	v[1] += rhs[1];
	v[2] += rhs[2];
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator-=(const nTYPE* rhs)
{
	v[0] -= rhs[0];
	v[1] -= rhs[1];
	v[2] -= rhs[2];
	return *this;
}

/* with scalars */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator+=(const nTYPE& value)
{
	v[0] += value;
	v[1] += value;
	v[2] += value;
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator-=(const nTYPE& value)
{
	v[0] -= value;
	v[1] -= value;
	v[2] -= value;
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator*=(const nTYPE& value)
{
	v[0] *= value;
	v[1] *= value;
	v[2] *= value;
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::operator/=(const nTYPE& value)
{
	v[0] /= value;
	v[1] /= value;
	v[2] /= value;
	return *this;
}

/* average */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::Average(const nTYPE* a, 
	const nTYPE* b)
{
	v[0] = 0.5*(a[0] + b[0]);
	v[1] = 0.5*(a[1] + b[1]);
	v[2] = 0.5*(a[2] + b[2]);
	return *this;
}

template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::Average(const nTYPE* a, 
	const nTYPE* b, const nTYPE* c)
{
	v[0] = (a[0] + b[0] + c[0])/3.0;
	v[1] = (a[1] + b[1] + c[1])/3.0;
	v[2] = (a[2] + b[2] + c[2])/3.0;
	return *this;
}


/* difference */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::Diff(const nTYPE* a, 
	const nTYPE* b)
{
	v[0] = a[0] - b[0];
	v[1] = a[1] - b[1];
	v[2] = a[2] - b[2];
	return *this;
}

/* vector cross product */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::Cross(const nTYPE* a, 
	const nTYPE* b)
{
	v[0] = a[1]*b[2] - a[2]*b[1];
	v[1] = a[2]*b[0] - a[0]*b[2];
	v[2] = a[0]*b[1] - a[1]*b[0];
	return *this;
}

/* fill the array with random numbers in the range [-1 1] */
template <class nTYPE>
inline void Vector3T<nTYPE>::Random(int seed)
{
	/* set random number seed */
	srand(seed);

	v[0] = nTYPE(rand() - RAND_MAX/2)/nTYPE(RAND_MAX);
	v[1] = nTYPE(rand() - RAND_MAX/2)/nTYPE(RAND_MAX);
	v[2] = nTYPE(rand() - RAND_MAX/2)/nTYPE(RAND_MAX);
}

/* inner product */
template <class nTYPE>
inline nTYPE Vector3T<nTYPE>::Dot(const nTYPE* a, const nTYPE* b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* L2 norm (between a and b) */
template <class nTYPE>
inline nTYPE Vector3T<nTYPE>::Norm(void) const
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

template <class nTYPE>
inline nTYPE Vector3T<nTYPE>::Norm(const nTYPE* a, const nTYPE* b)
{
	nTYPE dx = a[0] - b[0];
	nTYPE dy = a[1] - b[1];
	nTYPE dz = a[2] - b[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

/* linear combinations <- a*A + b*B */
template <class nTYPE>
inline Vector3T<nTYPE>& Vector3T<nTYPE>::
	Combine(const nTYPE& a, const nTYPE* A,
	        const nTYPE& b, const nTYPE* B)
{
	v[0] = a*A[0] + b*B[0];
	v[1] = a*A[1] + b*B[1];
	v[2] = a*A[2] + b*B[2];
	return *this;
}

} // namespace Tahoe 
#endif /* _VECTOR_3_T_H_ */
