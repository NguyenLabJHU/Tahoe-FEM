/* $Id: nArrayT.h,v 1.24 2004-03-16 05:37:14 paklein Exp $ */
/* created: paklein (05/23/1997) */
#ifndef _NARRAY_T_H_
#define _NARRAY_T_H_

/* ANSI headers */
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

#include "toolboxConstants.h"

/* base class */
#include "ArrayT.h"

namespace Tahoe {

/* forward declarations */
template <class nTYPE> class OutputProxyT;

/** templated base class for arrays of number types. The number type must
 * have the following mathematical operators defined:\n
 * (1) {+=, -=, *=, /=} : MUST return references to this\n
 * (2) { =, >>, < , > , fabs}\n
 * Additionally, the number type must allow assignment to 0.
 * The class extends the basic memory operations inherited from ArrayT
 * without adding any virtual functions. */
template <class nTYPE>
class nArrayT: public ArrayT<nTYPE>
{
public:

	/** \name constructors */
	/*@{*/
	/** default constructor. Constructs a zero length array. */
	nArrayT(void);

	/** construct an array of the specified length. The values
	 * in the array are not initialized.
	 * \param length length of dynamically allocated space */
	explicit nArrayT(int length);

	/** construct an alias
	 * \param length logical size of the array
	 * \param TYPEPtr pointer to the memory to use */
	nArrayT(int length, const nTYPE* TYPEPtr);

	/** copy constructor */
	nArrayT(const nArrayT& source);
	/*@}*/

	/** \name assignment operators */
	/*@{*/
	nArrayT<nTYPE>& operator=(const nArrayT& RHS); /**< assignment operator. Redimensions the array too match the source. */
	nArrayT<nTYPE>& operator=(const nTYPE* pRHS);  /**< assignment operator. Copy as many values as fit. */
	nArrayT<nTYPE>& operator=(const nTYPE& value); /**< set all elements in the array to value */
	/*@}*/

	/** \name addition operators */
	/*@{*/
	nArrayT<nTYPE>& operator+=(const nArrayT& RHS); /**< element-by-element addition with RHS */
	nArrayT<nTYPE>& operator+=(const nTYPE* pRHS);  /**< element-by-element addition with pRHS (without range checking). */
	nArrayT<nTYPE>& operator+=(const nTYPE& value); /**< add value to all elements */
	/*@}*/

	/** \name subtraction operators */
	/*@{*/
	nArrayT<nTYPE>& operator-=(const nArrayT& RHS); /**< element-by-element subtraction with RHS */
	nArrayT<nTYPE>& operator-=(const nTYPE* pRHS);  /**< element-by-element subtraction with pRHS (without range checking). */
	nArrayT<nTYPE>& operator-=(const nTYPE& value); /**< subtract value to all elements */
	/*@}*/

	/** \name multiplication operators */
	/*@{*/	
	nArrayT<nTYPE>& operator*=(const nArrayT& RHS); /**< element-by-element multiplication by RHS */
	nArrayT<nTYPE>& operator*=(const nTYPE& value); /**< multiply all elements by value */
	/*@}*/
	
	/** \name division operators */
	/*@{*/
	nArrayT<nTYPE>& operator/=(const nArrayT& RHS); /**< element-by-element division by RHS */		  	
	nArrayT<nTYPE>& operator/=(const nTYPE& value); /**< divide all elements by value */
	/*@}*/

	/** (post-)increment all elements in the array */
	nArrayT<nTYPE>& operator++(int);

	/** (post-)decrement all elements in the array */
	nArrayT<nTYPE>& operator--(int);
	
	/** return the sum of all elements in the array */
	nTYPE Sum(void) const;

	/** return the abs sum of all elements in the array */
	nTYPE AbsSum(void) const;

	/** return the average of the elements in the array */
	nTYPE Average(void) const;

	/** return the product of the elements in the array */
	nTYPE Product(void) const;

	/** \name inner products */
	/*@{*/
	static nTYPE Dot(const nArrayT<nTYPE>& A1, const nArrayT<nTYPE>& A2);
	static nTYPE Dot(const nTYPE* A1, const nArrayT<nTYPE>& A2);
	static nTYPE Dot(const nArrayT<nTYPE>& A1, const nTYPE* A2);
	static nTYPE Dot(const nTYPE* A1, const nTYPE* A2, int length);
	/*@}*/
	
	/** norm of the difference of two nArrayT's */
	static nTYPE Distance(const nArrayT<nTYPE>& A1, const nArrayT<nTYPE>& A2);

	/** \name max and min functions */
	/*@{*/
	nTYPE Max(void) const;          /**< return the maximum value in the array */
	nTYPE Max(int& position) const; /**< return the maximum value in the array and its position */
	nTYPE Min(void) const;          /**< return the minimum value in the array */
	nTYPE Min(int& position) const; /**< return the minimum value in the array and its position */
	nTYPE AbsMax(void) const;       /**< return the value with maximum absolute value */
	nTYPE AbsMin(void) const;       /**< return the value with minimum absolute value */
	void MinMax(nTYPE& min, nTYPE& max, bool positive_only = false) const; /**< return min and max values */
	void AbsMinMax(nTYPE& absmin, nTYPE& absmax) const; /**< return values values with the min and max absolute value */
	/*@}*/
	
	/** set all values with an absolute value smaller than tolerance to 0.0 */
	void Chop(double tolerance = kSmall);

	/** set array value to its offset from the beginning of the array, type cast
	 * the number type of the array */
	void SetValueToPosition(void);

	/** \name sorting */
	/*@{*/
	void SortAscending(void);
	void SortAscending(ArrayT<int>& map);
	void SortDescending(void); // not efficient
	/*@}*/

	/* commonly used operations:
	 *
	 *      this  = scale*RHS;
	 *		this += scale*RHS;
	 */		
	void SetToScaled(const nTYPE& scale, const nArrayT& RHS);
	void SetToScaled(const nTYPE& scale, const nTYPE* RHS);
	void AddScaled(const nTYPE& scale, const nArrayT& RHS);
	void AddScaled(const nTYPE& scale, const nTYPE* RHS);

	/** sum and difference of nArrayT's */
	void SumOf(const nArrayT& A, const nArrayT& B);
	void DiffOf(const nArrayT& A, const nArrayT& B);

	/** sum and difference arrays. Arrays assumed to be the same
	 * length as this. */
	void SumOf(const nTYPE* a, const nTYPE* b);
	void DiffOf(const nTYPE* a, const nTYPE* b);
	
	/* linear combinations:
	 *
	 *		this  = a*A + b*B
	 *		this  = a*A + b*B + c*C
	 *		this += a*A + b*B
	 *      this += a_1*A_1 + a_2*A_2 + ... + a_n-1*A_n-1
	 */
	void SetToCombination(const nTYPE& a, const nArrayT& A,
	                      const nTYPE& b, const nArrayT& B);
	void SetToCombination(const nTYPE& a, const nArrayT& A,
	                      const nTYPE& b, const nArrayT& B,
	                      const nTYPE& c, const nArrayT& C);
	void AddCombination(const nTYPE& a, const nArrayT& A,
	                    const nTYPE& b, const nArrayT& B);
	void AddCombination(const nArrayT& a, const ArrayT<nArrayT<nTYPE>*>& A);
	void AddCombination(const nTYPE& a, const nArrayT& A);

	/** fill the array with random numbers in the range [-1 1] */
	void Random(int seed = 1);
	
	/** \name I/O methods */
	/*@{*/
	void WriteWithFormat(ostream& out, int width, int prec,
		int wrapat, int tab = 0) const;

	/* write (with line wrapping) */
	void WriteNoWrap(ostream& out) const;	
	void WriteNoWrapTight(ostream& out) const;
	void WriteWrapped(ostream& out, int line_count, int tab = 0) const;
	void WriteWrappedTight(ostream& out, int line_count) const;

	/* proxies for "<<" */
	OutputProxyT<nTYPE> no_wrap(void) const {
		return OutputProxyT<nTYPE>(OutputProxyT<nTYPE>::kNoWrap, *this, 0, 0); };

	OutputProxyT<nTYPE> no_wrap_tight(void) const {
		return OutputProxyT<nTYPE>(OutputProxyT<nTYPE>::kNoWrapTight, *this, 0, 0); };

	OutputProxyT<nTYPE> wrap(int line_count, int tab = 0) const {
		return OutputProxyT<nTYPE>(OutputProxyT<nTYPE>::kWrap, *this, line_count, tab); };

	OutputProxyT<nTYPE> wrap_tight(int line_count, int tab = 0) const {
		return OutputProxyT<nTYPE>(OutputProxyT<nTYPE>::kWrapTight, *this, line_count, tab); };
	/*@}*/
};


/* I/O operators */
template <class nTYPE>
istream& operator>>(istream& in, nArrayT<nTYPE>& array)
{
	nTYPE* p = array.Pointer();
	for (int i = 0; i < array.Length(); i++)
		in >> *p++;
	return in;
}

template <class nTYPE>
ostream& operator<<(ostream& out, const nArrayT<nTYPE>& array)
{
	const nTYPE* p = array.Pointer();
	int width = OutputWidth(out, p);
	for (int i = 0; i < array.Length(); i++)
	{
		if (i > 0) out << '\n';
		out << setw(width) << *p++;
	}
	return out;
}

/* output formatters proxy for use with << */
template <class TYPE>
class OutputProxyT
{
public:

	/* formats */
	enum FormatT {kNoWrap, kNoWrapTight, kWrap, kWrapTight};

	/* constructor */
	OutputProxyT(FormatT format, const nArrayT<TYPE>& array, int line_count, int tab):
		format_(format),
		array_(array),
		line_count_(line_count),
		tab_(tab) {};

public:

	FormatT format_;
	const nArrayT<TYPE>& array_;
	int line_count_;
	int tab_;
};

/* output operator */
template <class TYPE>
ostream& operator<<(ostream& out, const OutputProxyT<TYPE>& proxy)
{
	switch (proxy.format_)
	{
		case OutputProxyT<TYPE>::kNoWrap:
			proxy.array_.WriteNoWrap(out);			
			break;
		case OutputProxyT<TYPE>::kNoWrapTight:
			proxy.array_.WriteNoWrapTight(out);			
			break;
		case OutputProxyT<TYPE>::kWrap:
			proxy.array_.WriteWrapped(out, proxy.line_count_, proxy.tab_);
			break;
		case OutputProxyT<TYPE>::kWrapTight:
			proxy.array_.WriteWrappedTight(out, proxy.line_count_);
			break;
		default:
			ExceptionT::GeneralFail("operator<<OutputProxyT", "unknown format: %d", proxy.format_);
	}
	return out;
}

/* output formatters - for int's and double's */
inline int OutputWidth(ostream& out, const int* junk)
{
#pragma unused(junk)
#pragma unused(out)
	return kIntWidth;
};

inline int OutputWidth(ostream& out, const float* junk)
{
#pragma unused(junk)
	return out.precision() + kDoubleExtra;
};

inline int OutputWidth(ostream& out, const double* junk)
{
#pragma unused(junk)
	return out.precision() + kDoubleExtra;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(void) { }

template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(int length): ArrayT<nTYPE>(length) { }

template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(int length, const nTYPE* MATHTYPEPtr):
	ArrayT<nTYPE>(length, MATHTYPEPtr) { }

template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(const nArrayT& source):
	ArrayT<nTYPE>(source) { }

/* write with line wrapping */
template <class nTYPE>
void nArrayT<nTYPE>::WriteNoWrap(ostream& out) const
{
	const nTYPE* p = Pointer();
	int width = OutputWidth(out, p);
	for (int i = 0; i < Length(); i++)
		out << setw(width) << p[i];
}

template <class nTYPE>
void nArrayT<nTYPE>::WriteNoWrapTight(ostream& out) const
{
	const nTYPE* p = Pointer();
	for (int i = 0; i < Length(); i++)
		out << " " << p[i];
}

template <class nTYPE>
void nArrayT<nTYPE>::WriteWrapped(ostream& out, int linecount, int tab) const
{
	const nTYPE*  p = Pointer();
	int width = OutputWidth(out, p);
	int count = 0;
	for (int i = 0; i < Length(); i++)
	{
		/* wrap */
		if (count == linecount)
		{
			out << '\n';
			count = 0;

			/* tab out */
			if (tab > 0) out << setw(tab) << " ";
		}
		count++;
		out << setw(width) << *p++;		
	}
}

/* print with line wrapping */
template <class nTYPE>
void nArrayT<nTYPE>::WriteWrappedTight(ostream& out, int linecount) const
{
	const nTYPE*  p = Pointer();
//	int width = OutputWidth(out, p);
	int count = 0;
	for (int i = 0; i < Length(); i++)
	{
		/* wrap */
		if (count == linecount)
		{
			out << '\n';
			count = 0;
		}
		else if (count > 0)
			out << " ";
		count++;

		out << *p++;		
	}
}

/* copy/assignment operators - by a scalar or element by element */
template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator=(const nArrayT& RHS)
{
	/* inherited */
	ArrayT<nTYPE>::operator=(RHS);
	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator=(const nTYPE* pRHS)
{
#if __option (extended_errorcheck)
	/* check */
	if (!pRHS) ExceptionT::GeneralFail("nArrayT<nTYPE>::operator=", "pointer is NULL");
#endif

	/* no copies to self */
	if (pRHS != Pointer()) MemCopy(Pointer(), pRHS, Length());

	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator=(const nTYPE& value)
{
	/* inherited */
	ArrayT<nTYPE>::operator=(value);
	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator+=(const nArrayT& RHS)
{
#if __option (extended_errorcheck)
	/* dimension checks */
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch("nArrayT<nTYPE>::operator+=");
#endif
	return operator+=(RHS.Pointer());
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator+=(const nTYPE* pRHS)
{
	nTYPE* pthis = Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ += *pRHS++;
		
	return *this ;
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator+=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ += value;
		
	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator-=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch("nArrayT<nTYPE>::operator-=");
#endif
	return operator-=(RHS.Pointer());
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator-=(const nTYPE* pRHS)
{
	nTYPE* pthis = Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ -= *pRHS++;
		
	return *this;
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator-=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ -= value;
		
	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator*=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch("nArrayT<nTYPE>::operator*=");
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pRHS  = RHS.Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ *= *pRHS++;
		
	return *this;
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator*=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ *= value;
		
	return *this;	
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator/=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch("nArrayT<nTYPE>::operator/=");
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pRHS  = RHS.Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ /= *pRHS++;
		
	return *this;
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator/=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ /= value;
		
	return *this;	
}

/* increment/decriment all */
template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator++(int)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
	{
		(*pA)++;
		pA++;
	}
	return *this;
}

template <class nTYPE>
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator--(int)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
	{
		(*pA)--;
		pA++;
	}
	return *this;
}

/* sum, average, and product */
template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Sum(void) const
{
	register nTYPE sum = nTYPE(0.0);
	const nTYPE* p = Pointer();
	for (int i = 0; i < fLength; i++)
		sum += *p++;
	return sum;
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::AbsSum(void) const
{
	register nTYPE sum = nTYPE(0.0);
	const nTYPE* p = Pointer();
	for (int i = 0; i < fLength; i++)
		sum += fabs(*p++);
	return sum;
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Average(void) const
{
	return Sum()/fLength;
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Product(void) const
{
	register nTYPE product = nTYPE(1.0);
	const nTYPE* p = Pointer();
	for (int i = 0; i < fLength; i++)
		product *= *p++;
	return product;
}

/* max and min */
template <class nTYPE>
nTYPE nArrayT<nTYPE>::Max(void) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE max = *pthis++;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis > max) max = *pthis;
		pthis++;
	}

	return max;
}

template <class nTYPE>
nTYPE nArrayT<nTYPE>::Max(int& position) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE max = *pthis++;
	position = 0;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis > max)
		{
			max = *pthis;
			position = i;
		}
		pthis++;
	}

	return max;
}

template <class nTYPE>
nTYPE nArrayT<nTYPE>::Min(void) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE min = *pthis++;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis < min) min = *pthis;
		pthis++;
	}

	return min;
}

template <class nTYPE>
nTYPE nArrayT<nTYPE>::Min(int& position) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE  min   = *pthis++;
	position = 0;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis < min)
		{
			min = *pthis;
			position = i;
		}
		pthis++;
	}

	return min;
}

template <class nTYPE>
nTYPE nArrayT<nTYPE>::AbsMax(void) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE abs, max = fabs(*pthis++);
	for (int i = 1; i < Length(); i++)
	{
		abs = fabs(*pthis++);
		if (abs > max) max = abs;
	}

	return max;
}

template <class nTYPE>
nTYPE nArrayT<nTYPE>::AbsMin(void) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	const nTYPE* pthis = Pointer();
	nTYPE abs, min = fabs(*pthis++);
	for (int i = 1; i < Length(); i++)
	{
		abs = fabs(*pthis++);
		if (abs < min) min = abs;
	}

	return min;
}

template <class nTYPE>
void nArrayT<nTYPE>::MinMax(nTYPE& min, nTYPE& max,
	bool positive_only) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	/* ignore negative numbers */
	if (positive_only)
	{
		const nTYPE* pthis = Pointer();
		max = 0;
		min = *pthis++;
		for (int i = 1; i < Length(); i++)
		{
			/* try to keep positive */
			if (min < 0) min = *pthis;
		
			/* same */
			if (*pthis < min && *pthis >= 0)
				min = *pthis;
			else if (*pthis > max)
				max = *pthis;
				
			pthis++;
		}
		
		if (min < 0) min = 0;
	}
	else
	{
		const nTYPE* pthis = Pointer();
		min = *pthis++;
		max = min;
		for (int i = 1; i < Length(); i++)
		{
			if (*pthis < min)
				min = *pthis;
			else if (*pthis > max)
				max = *pthis;
				
			pthis++;
		}
	}
}

template <class nTYPE>
void nArrayT<nTYPE>::AbsMinMax(nTYPE& absmin, nTYPE& absmax) const
{
#if __option(extended_errorcheck)
	if (!fArray) ExceptionT::GeneralFail();
#endif

	nTYPE abs;
	const nTYPE* pthis = Pointer();
	absmax = absmin = fabs(*pthis++);
	for (int i = 1; i < Length(); i++)
	{
		abs = fabs(*pthis++);

		if (abs < absmin)
			absmin = abs;
		else if (abs > absmin)
			absmax = abs;
	}
}

/* removing small values */
template <class nTYPE>
void nArrayT<nTYPE>::Chop(double tolerance)
{
	nTYPE* pthis = Pointer();
	for (int i = 0; i < fLength; i++)
	{
		if (fabs(*pthis) < tolerance)
			*pthis = 0;
	
		pthis++;	
	}
}

/* set array value to its position in the array */
template <class nTYPE>
void nArrayT<nTYPE>::SetValueToPosition(void)
{
	nTYPE* p = Pointer();
	for (int i = 0; i < Length(); i++)
		*p++ = nTYPE(i);
}

/* sorting */
template <class nTYPE>
void nArrayT<nTYPE>::SortAscending(void)
{
	nTYPE* list = Pointer();
	int N = Length();

	int l, r, j, i, flag;
	nTYPE RR, K;

	if (N <= 1) return;

	l   = N / 2 + 1;
	r   = N - 1;
	l   = l - 1;
	RR  = list[l - 1];
	K   = list[l - 1];

	while (r != 0) {
		j = l;
		flag = 1;

		while (flag == 1) {
	i = j;
	j = j + j;

			if (j > r + 1)
				flag = 0;
			else {
				if (j < r + 1)
					if (list[j] > list[j - 1]) j = j + 1;

				if (list[j - 1] > K) {
					list[i - 1] = list[j - 1];
				}
				else {
					flag = 0;
				}
			}
		}

		list[ i - 1] = RR;

		if (l == 1) {
			RR  = list [r];

			K = list[r];
			list[r ] = list[0];
			r = r - 1;
		}
		else {
			l   = l - 1;
			RR  = list[l - 1];
			K   = list[l - 1];
		}
	}

	list[0] = RR;
}

/* sort this and map, by values in map */
template <class nTYPE>
inline void nArrayT<nTYPE>::SortAscending(ArrayT<int>& map)
{
	Tahoe::SortAscending(map, *this);
}

/* AZTEC az_sort.c: sort both arrays by the values in the master array */
template <class masterTYPE, class slaveTYPE>
static void SortAscending(ArrayT<masterTYPE>& master, ArrayT<slaveTYPE>& slave)
{
	/* check */
	if (master.Length() < slave.Length()) ExceptionT::SizeMismatch();
	int N = master.Length();

	/* local variables */
	int    l, r, j, i, flag;
	masterTYPE RR, K;
	slaveTYPE  RR2;

	masterTYPE* list = master.Pointer();
	slaveTYPE* list2 = slave.Pointer();

if (N <= 1) return;

l   = N / 2 + 1;
r   = N - 1;
l   = l - 1;
RR  = list[l - 1];
K   = list[l - 1];

RR2 = list2[l - 1];
while (r != 0) {
j = l;
flag = 1;

while (flag == 1) {
i = j;
j = j + j;

if (j > r + 1)
flag = 0;
else {
if (j < r + 1)
if (list[j] > list[j - 1]) j = j + 1;

if (list[j - 1] > K) {
list[ i - 1] = list[ j - 1];
list2[i - 1] = list2[j - 1];
}
else {
flag = 0;
}
}
}

list[ i - 1] = RR;
list2[i - 1] = RR2;

if (l == 1) {
RR  = list [r];
RR2 = list2[r];

K = list[r];
list[r ] = list[0];
list2[r] = list2[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR2 = list2[l - 1];
K   = list[l - 1];
}
}

list[ 0] = RR;
list2[0] = RR2;
} /* AZ_sort */

template <class nTYPE>
void nArrayT<nTYPE>::SortDescending(void)
{
	nTYPE temp;

	int swaps = 1;
	while (swaps > 0)
	{
		swaps = 0;
		
		nTYPE* p = Pointer();
		for (int i = 1; i < fLength; i++)
		{
			if (*p < *(p+1))
			{
				swaps++;
			
				temp = *p;
				*p     = *(p+1);
				*(p+1) = temp;
			}
			
			p++;
		}	
	}
}

/* compute the inner product of a1 and a2 */
template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Dot(const nTYPE* A1, const nTYPE* A2, int length)
{
	register nTYPE dot = 0.0;
	register nTYPE temp;
	for (int i = 0; i < length; i++) {
		temp  = (*A1++);
		temp *= (*A2++);
		dot += temp;
	}
	return dot;
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Dot(const nArrayT<nTYPE>& A1, const nArrayT<nTYPE>& A2)
{
/* dimension check */
#if __option (extended_errorcheck)
	if (A1.Length() != A2.Length()) ExceptionT::SizeMismatch("nArrayT<nTYPE>::Dot");
#endif
	return Dot(A1.Pointer(), A2.Pointer(), A1.Length());
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Dot(const nTYPE* A1, const nArrayT<nTYPE>& A2)
{
	return Dot(A1, A2.Pointer(), A2.Length());
}

template <class nTYPE>
inline nTYPE nArrayT<nTYPE>::Dot(const nArrayT<nTYPE>& A1, const nTYPE* A2)
{
	return Dot(A1.Pointer(), A2, A1.Length());
}


/* distance */
template <class nTYPE>
nTYPE nArrayT<nTYPE>::Distance(const nArrayT<nTYPE>& A1, 
	const nArrayT<nTYPE>& A2)
{
/* dimension check */
#if __option (extended_errorcheck)
	if (A1.Length() != A2.Length()) ExceptionT::SizeMismatch();
#endif

	const nTYPE* p1 = A1.Pointer();
	const nTYPE* p2 = A2.Pointer();
	
	register nTYPE dot = 0.0;
	register nTYPE temp;
	
	int length = A1.Length();
	for (int i = 0; i < length; i++)
	{
		temp  = (*p1++);
		temp -= (*p2++);
		temp *= temp;
		dot += temp;
	}
	return sqrt(dot);
}

/* commonly used operations;
*
*      this =  scale*RHS;
*		this += scale*RHS;
*
*  Note: returns a reference to (*this) for chaining operations */		
template <class nTYPE>
inline void nArrayT<nTYPE>::SetToScaled(const nTYPE& scale, const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch();
#endif
	SetToScaled(scale, RHS.Pointer());
}

template <class nTYPE>
void nArrayT<nTYPE>::SetToScaled(const nTYPE& scale, const nTYPE* RHS)
{
	nTYPE* pthis = Pointer();
	nTYPE temp;
	for (int i = 0; i < fLength; i++)
	{
		temp  = scale;
		temp *= *RHS++;	//multi-step needed incase pthis == RHS
		*pthis++ = temp;
	}
}

template <class nTYPE>
inline void nArrayT<nTYPE>::AddScaled(const nTYPE& scale, const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) ExceptionT::SizeMismatch();
#endif

	AddScaled(scale, RHS.Pointer());
}

template <class nTYPE>
void nArrayT<nTYPE>::AddScaled(const nTYPE& scale, const nTYPE* RHS)
{
	nTYPE* pthis = Pointer();
	nTYPE temp;
	for (int i = 0; i < fLength; i++)
	{
		temp  = scale;
		temp *= *RHS++;
		*pthis++ += temp;
	}
}

/* sum (A + B) and difference (A - B) of nArrayT's */
template <class nTYPE>
inline void nArrayT<nTYPE>::SumOf(const nArrayT& A, const nArrayT& B)
{
#if __option (extended_errorcheck)
	/* dimension checks */
	if (fLength != A.fLength || fLength != B.fLength) ExceptionT::SizeMismatch();
#endif

	/* call pointer version */
	SumOf(A.Pointer(), B.Pointer());
}

template <class nTYPE>
inline void nArrayT<nTYPE>::SumOf(const nTYPE* pA, const nTYPE* pB)
{
	nTYPE* pthis = Pointer();
	for (int i = 0; i < fLength; i++)
	{
		*pthis  = *pA++;
		*pthis += *pB++;
		pthis++;
	}
}
	
template <class nTYPE>
inline void nArrayT<nTYPE>::DiffOf(const nArrayT& A, const nArrayT& B)
{
#if __option (extended_errorcheck)
	/* dimension checks */
	if (fLength != A.fLength || fLength != B.fLength) ExceptionT::SizeMismatch();
#endif

	/* call pointer version */
	DiffOf(A.Pointer(), B.Pointer());
}

template <class nTYPE>
void nArrayT<nTYPE>::DiffOf(const nTYPE* pA, const nTYPE* pB)
{
	nTYPE* pthis = Pointer();
	for (int i = 0; i < fLength; i++)
	{
		*pthis  = *pA++;
		*pthis -= *pB++;
		pthis++;
	}
}

/* linear combinations */
template <class nTYPE>
void nArrayT<nTYPE>::SetToCombination(const nTYPE& a, const nArrayT& A,
	const nTYPE& b, const nArrayT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength || fLength != B.fLength) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pA = A.Pointer();
	const nTYPE* pB = B.Pointer();

	register nTYPE temp;
	
	for (int i = 0; i < fLength; i++)
	{
		temp  = a;
		temp *= *pA++;
		*pthis = temp;
	
		temp  = b;
		temp *= *pB++;
		*pthis++ += temp;
	}
}

template <class nTYPE>
void nArrayT<nTYPE>::SetToCombination(
	const nTYPE& a, const nArrayT& A,
	const nTYPE& b, const nArrayT& B,
	const nTYPE& c, const nArrayT& C)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength ||
	    fLength != B.fLength ||
	    fLength != C.fLength) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pA = A.Pointer();
	const nTYPE* pB = B.Pointer();
	const nTYPE* pC = C.Pointer();

	register nTYPE temp;
	
	for (int i = 0; i < fLength; i++)
	{
		temp  = a;
		temp *= *pA++;
		*pthis = temp;
	
		temp  = b;
		temp *= *pB++;
		*pthis += temp;

		temp  = c;
		temp *= *pC++;
		*pthis++ += temp;
	}
}

template <class nTYPE>
void nArrayT<nTYPE>::AddCombination(const nTYPE& a, const nArrayT& A,
	const nTYPE& b, const nArrayT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength || fLength != B.fLength) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pA = A.Pointer();
	const nTYPE* pB = B.Pointer();

	register nTYPE temp;
	
	for (int i = 0; i < fLength; i++)
	{
		temp  = a;
		temp *= *pA++;
		*pthis += temp;
	
		temp  = b;
		temp *= *pB++;
		*pthis++ += temp;
	}
}

template <class nTYPE>
void nArrayT<nTYPE>::AddCombination(const nTYPE& a, const nArrayT& A)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength) ExceptionT::SizeMismatch();
#endif

	nTYPE* pthis = Pointer();
	const nTYPE* pA = A.Pointer();

	register nTYPE temp;
	
	for (int i = 0; i < fLength; i++)
	{
		temp  = a;
		temp *= *pA++;
		*pthis += temp;
	}
}

template <class nTYPE>
void nArrayT<nTYPE>::AddCombination(const nArrayT& a,
	const ArrayT<nArrayT<nTYPE>*>& A)
{
/* dimension of sum */
#if __option (extended_errorcheck)
	if (a.Length() != A.Length()) ExceptionT::SizeMismatch();
#endif

	/* sum */
	for (int i = 0; i < a.Length(); i++)
	{
		/* check each term in sum */
		if (fLength != A[i]->Length()) ExceptionT::SizeMismatch();
	
		nTYPE* pthis = Pointer();
		nTYPE* pA    = A[i]->Pointer();
	
		register nTYPE temp;
		
		for (int j = 0; j < fLength; j++)
		{
			temp  = a[i];
			temp *= *pA++;
			*pthis++ += temp;
		}
	}
}

/* fill the array with random numbers in the range [-1 1] */
template <class nTYPE>
void nArrayT<nTYPE>::Random(int seed)
{
	/* set random number seed */
	srand(seed);

	nTYPE* p = Pointer();
	for (int i = 0; i < Length(); i++)
		*p++ = nTYPE(rand() - RAND_MAX/2)/nTYPE(RAND_MAX);
}

/* output */
template <class nTYPE>
void nArrayT<nTYPE>::WriteWithFormat(ostream& out, int width, int prec,
	int wrapat, int tab) const
{
	int currprec = out.precision();
	out.precision(prec);
	
	int     count = 0;
	const nTYPE* pthis = Pointer();
	
	for (int i = 0; i < Length(); i++)
	{
		out << setw(width) << *pthis++;
		
		/* wrap */
		if (++count == wrapat)
		{
			out << '\n';
			count = 0;
			
			/* tab out */
			if (tab > 0) out << setw(tab) << " ";
		}
	}

	if (count != 0) out << '\n';

	out.precision(currprec);
}

} /* namespace Tahoe */

#endif /* _NARRAY_T_H_ */
