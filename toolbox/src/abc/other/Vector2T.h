/* $Id: Vector2T.h,v 1.4 2005-03-12 08:36:10 paklein Exp $ */
#ifndef _VECTOR_2_T_H_
#define _VECTOR_2_T_H_

/* Environmental */
#include "Environment.h"

namespace Tahoe {

/** utility class for 2D vector functions
 * NOTE: some functions do create temporary nTYPE's */
template <class nTYPE>
class Vector2T
{
  public:

	/* constructors */
	Vector2T(void);
	Vector2T(const nTYPE* source);

	/** \name accessors */
	/*@{*/
	nTYPE& operator[](int dex);
	const nTYPE& operator[](int dex) const;

	/* type conversion operators */
	operator const nTYPE*() const;
	operator nTYPE*();
	/*@}*/

	/* assignment operator */
	Vector2T<nTYPE>& operator=(const nTYPE* rhs);
	Vector2T<nTYPE>& operator=(const nTYPE& value);

	/* mathematical operators */
	Vector2T<nTYPE>& operator+=(const nTYPE* rhs);
	Vector2T<nTYPE>& operator-=(const nTYPE* rhs);

	/* with scalars */
	Vector2T<nTYPE>& operator+=(const nTYPE& value);
	Vector2T<nTYPE>& operator-=(const nTYPE& value);
	Vector2T<nTYPE>& operator*=(const nTYPE& value);
	Vector2T<nTYPE>& operator/=(const nTYPE& value);

	/* average */
	Vector2T<nTYPE>& Average(const nTYPE* a, const nTYPE* b);
	Vector2T<nTYPE>& Average(const nTYPE* a, const nTYPE* b, const nTYPE* c);
	
	/* difference */
	Vector2T<nTYPE>& Diff(const nTYPE* a, const nTYPE* b);

	/* vector cross product */
	Vector2T<nTYPE>& Cross(const nTYPE* a);

	/* inner product */
	static nTYPE Dot(const nTYPE* a, const nTYPE* b);

	/* L2 norm (between a and b) */
	nTYPE Norm(void) const;
	static nTYPE Norm(const nTYPE* a, const nTYPE* b);

	/* linear combinations <- a*A + b*B */
	Vector2T<nTYPE>& Combine(const nTYPE& a, const nTYPE* A,
	                         const nTYPE& b, const nTYPE* B);

private:

	nTYPE v[2];
};

/* inlines */
template <class nTYPE>
inline Vector2T<nTYPE>::Vector2T(void)
{ 
	v[0] = 0;
	v[1] = 0;
}

template <class nTYPE>
inline Vector2T<nTYPE>::Vector2T(const nTYPE* source)
{
	operator=(source);
}

/* accessor */
template <class nTYPE>
inline nTYPE& Vector2T<nTYPE>::operator[](int dex) const
{
#if __option(extended_errorcheck)
	if (dex < 0 || dex > 2) throw eOutOfRange;
#endif
	return v[dex];
}

/* type conversion operators */
template <class nTYPE>
inline Vector2T<nTYPE>::operator const nTYPE*() const { return (const nTYPE*) v; }

template <class nTYPE>
inline Vector2T<nTYPE>::operator nTYPE*() { return (nTYPE*) v; }

/* assignment operator */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator=(const nTYPE* rhs)
{
	v[0] = rhs[0];
	v[1] = rhs[1];
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator=(const nTYPE& value)
{
	v[0] = value;
	v[1] = value;
	return *this;
}

/* mathematical operators */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator+=(const nTYPE* rhs)
{
	v[0] += rhs[0];
	v[1] += rhs[1];
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator-=(const nTYPE* rhs)
{
	v[0] -= rhs[0];
	v[1] -= rhs[1];
	return *this;
}

/* with scalars */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator+=(const nTYPE& value)
{
	v[0] += value;
	v[1] += value;
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator-=(const nTYPE& value)
{
	v[0] -= value;
	v[1] -= value;
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator*=(const nTYPE& value)
{
	v[0] *= value;
	v[1] *= value;
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::operator/=(const nTYPE& value)
{
	v[0] /= value;
	v[1] /= value;
	return *this;
}

/* average */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::Average(const nTYPE* a, 
	const nTYPE* b)
{
	v[0] = 0.5*(a[0] + b[0]);
	v[1] = 0.5*(a[1] + b[1]);
	return *this;
}

template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::Average(const nTYPE* a, 
	const nTYPE* b, const nTYPE* c)
{
	v[0] = (a[0] + b[0] + c[0])/3.0;
	v[1] = (a[1] + b[1] + c[1])/3.0;
	return *this;
}


/* difference */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::Diff(const nTYPE* a, 
	const nTYPE* b)
{
	v[0] = a[0] - b[0];
	v[1] = a[1] - b[1];
	return *this;
}

/* vector cross product */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::Cross(const nTYPE* a )
{ // this a X e3 
	v[0] = a[1];
	v[1] = - a[0];
	return *this;
}

/* inner product */
template <class nTYPE>
inline nTYPE Vector2T<nTYPE>::Dot(const nTYPE* a, const nTYPE* b)
{
	return a[0]*b[0] + a[1]*b[1] ;
}

/* L2 norm (between a and b) */
template <class nTYPE>
inline nTYPE Vector2T<nTYPE>::Norm(void) const
{
	return sqrt(v[0]*v[0] + v[1]*v[1] );
}

template <class nTYPE>
inline nTYPE Vector2T<nTYPE>::Norm(const nTYPE* a, const nTYPE* b)
{
	nTYPE dx = a[0] - b[0];
	nTYPE dy = a[1] - b[1];
	return sqrt(dx*dx + dy*dy);
}

/* linear combinations <- a*A + b*B */
template <class nTYPE>
inline Vector2T<nTYPE>& Vector2T<nTYPE>::
	Combine(const nTYPE& a, const nTYPE* A,
	        const nTYPE& b, const nTYPE* B)
{
	v[0] = a*A[0] + b*B[0];
	v[1] = a*A[1] + b*B[1];
	return *this;
}

} // namespace Tahoe 
#endif /* _VECTOR_2_T_H_ */
