/*
 * File: nArrayT.h - TEMPLATE
 *
 * Base class for arrays of TYPE for which the following mathematical
 * operators have been defined: 
 *
 *   {+=, -=, *=, /=} : MUST return references to this
 *   { =, >>, < , > , fabs}
 *
 * And must allow assignment to 0.0. TYPE must also contain no virtual
 * functions, since the Clear() assumes all bytes can be set to 0.
 *
 */

/*
 * created      : PAK (05/23/97)
 * last modified: PAK (07/03/98)
 */

#ifndef _NARRAY_T_H_
#define _NARRAY_T_H_

/* ANSI headers */
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"

/* base class */
#include "ArrayT.h"

template <class nTYPE>
class nArrayT: public ArrayT<nTYPE>
{
  public:

	/* constructors */
	nArrayT(void);
	nArrayT(int length);
	nArrayT(int length, nTYPE* TYPEPtr);
	nArrayT(const nArrayT& source); 

	/* write with line wrapping */
	void WriteWrapped(ostream& out, int linecount, int tab = 0) const; 
	void WriteWrappedTight(ostream& out, int linecount) const; 

  	/* copy/assignment operators - by a scalar or element by element */  	 
  	nArrayT<nTYPE>& operator=(const nArrayT& RHS);
  	nArrayT<nTYPE>& operator=(const nTYPE& value);

  	nArrayT<nTYPE>& operator+=(const nArrayT& RHS);
  	nArrayT<nTYPE>& operator+=(const nTYPE& value);

  	nArrayT<nTYPE>& operator-=(const nArrayT& RHS);
	nArrayT<nTYPE>& operator-=(const nTYPE& value);

  	nArrayT<nTYPE>& operator*=(const nArrayT& RHS);
  	nArrayT<nTYPE>& operator*=(const nTYPE& value);
  	
  	nArrayT<nTYPE>& operator/=(const nArrayT& RHS); 		  	 
  	nArrayT<nTYPE>& operator/=(const nTYPE& value);

	/* (post-)increment/decrement all */
  	nArrayT<nTYPE>& operator++(int);
  	nArrayT<nTYPE>& operator--(int);
 	
	/* sum, average, and product */
  	nTYPE Sum(void) const;
  	nTYPE Average(void) const;
  	nTYPE Product(void) const;

	/* max and min */
	nTYPE Max(void) const;
	nTYPE Max(int& position) const;
	nTYPE Min(void) const;
	nTYPE Min(int& position) const;
	nTYPE AbsMax(void) const;
	nTYPE AbsMin(void) const;
	void MinMax(nTYPE& min, nTYPE& max, bool positive_only = false) const;
	void AbsMinMax(nTYPE& absmin, nTYPE& absmax) const;
	
	/* removing small values */
	void Chop(double tolerance = kSmall);
	
	/* sorting */
	void SortAscending(void);
	void SortDescending(void);

	/* commonly used operations:
	 *
	 *      this  = scale*RHS;
	 *		this += scale*RHS;
	 */		
	void SetToScaled(const nTYPE& scale, const nArrayT& RHS);
	void AddScaled(const nTYPE& scale, const nArrayT& RHS);

	/* sum and difference of nArrayT's */
	void SumOf(const nArrayT& A, const nArrayT& B);
	void DiffOf(const nArrayT& A, const nArrayT& B);
	
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
	void AddCombination(const nArrayT& a, 
                        const ArrayT<nArrayT<nTYPE>*>& A);
                        
	/* fill the array with random numbers in the range [-1 1] */
	void Random(int seed = 1);
	
	/* output */
	void PrintWithFormat(ostream& out, int width, int prec,
		int wrapat, int tab = 0) const;
	
	/* I/O operators */
	friend istream& operator>>(istream& in, nArrayT<nTYPE>& array)
	{
		nTYPE* p = array.Pointer();
		for (int i = 0; i < array.Length(); i++)
			in >> *p++;
		return in;
	};

	friend ostream& operator<<(ostream& out, const nArrayT<nTYPE>& array)
	{
		nTYPE* p = array.Pointer();
		int width = OutputWidth(p, out);
		for (int i = 0; i < array.Length(); i++)
		{
			if (i > 0) out << '\n';
			out << setw(width) << *p++;
		}
		return out;
	};
};

/* output formatters - for int's and double's */
inline int OutputWidth(const int* junk, ostream& out)
{
#pragma unused(junk)
#pragma unused(out)
	return kIntWidth;
};

inline int OutputWidth(const float* junk, ostream& out)
{
#pragma unused(junk)
	return out.precision() + kDoubleExtra;
};

inline int OutputWidth(const double* junk, ostream& out)
{
#pragma unused(junk)
	return out.precision() + kDoubleExtra;
};

/*************************************************************************
 *
 * Implementation
 *
 *************************************************************************/

/* inner product */
template <class nTYPE>
nTYPE Dot(const nArrayT<nTYPE>& A1, 
	                    const nArrayT<nTYPE>& A2);

/* constructor */
template <class nTYPE> 
inline nArrayT<nTYPE>::nArrayT(void) { }

template <class nTYPE> 
inline nArrayT<nTYPE>::nArrayT(int length): ArrayT<nTYPE>(length) { }

template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(int length, nTYPE* MATHTYPEPtr):
	ArrayT<nTYPE>(length, MATHTYPEPtr) { }

template <class nTYPE>
inline nArrayT<nTYPE>::nArrayT(const nArrayT& source):
	ArrayT<nTYPE>(source) { }

/* print with line wrapping */
template <class nTYPE>
void nArrayT<nTYPE>::WriteWrapped(ostream& out, int linecount, int tab) const
{
	nTYPE*  p = Pointer();
	int width = OutputWidth(p, out);
	int count = 0;
	
	for (int i = 0; i < Length(); i++)
	{
		out << setw(width) << *p++;
		
		/* wrap */
		if (++count == linecount)
		{
			out << '\n';
			count = 0;

			/* tab out */
			if (tab > 0 && i < Length() - 1)
				out << setw(tab) << " ";
		}
	}

	if (count != 0) out << '\n';
}

/* print with line wrapping */
template <class nTYPE>
void nArrayT<nTYPE>::WriteWrappedTight(ostream& out, int linecount) const
{
	nTYPE*  p = Pointer();
	int width = OutputWidth(p, out);
	int count = 0;
	for (int i = 0; i < Length(); i++)
	{
		out << *p++;
		
		/* wrap */
		if (++count == linecount)
		{
			out << '\n';
			count = 0;
		}
		else
			out << " ";
	}

	if (count != 0) out << '\n';
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
inline nArrayT<nTYPE>& nArrayT<nTYPE>::operator=(const nTYPE& value)
{
	/* inherited */
	ArrayT<nTYPE>::operator=(value);
	return *this;	
}

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator+=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ += *pRHS++;
		
	return *this ;
} 

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator+=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ += value;
		
	return *this;	
}

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator-=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();

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
nArrayT<nTYPE>& nArrayT<nTYPE>::operator*=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ *= *pRHS++;
		
	return *this;
} 

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator*=(const nTYPE& value)
{
	nTYPE* pA = Pointer();
	for (int i = 0; i < fLength; i++)
		*pA++ *= value;
		
	return *this;	
}

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator/=(const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();
	for (int i = 0; i < fLength; i++)
		*pthis++ /= *pRHS++;
		
	return *this;
} 

template <class nTYPE> 
nArrayT<nTYPE>& nArrayT<nTYPE>::operator/=(const nTYPE& value)
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
nTYPE nArrayT<nTYPE>::Sum(void) const
{
	register nTYPE sum = nTYPE(0.0);
	nTYPE* p = Pointer();
	
	for (int i = 0; i < fLength; i++)
		sum += *p++;

	return sum;
}
 
template <class nTYPE> 
inline nTYPE nArrayT<nTYPE>::Average(void) const
{
	return Sum()/fLength;
}

template <class nTYPE> 
nTYPE nArrayT<nTYPE>::Product(void) const
{
	register nTYPE product = nTYPE(1.0);
	nTYPE* p = Pointer();
	for (int i = 0; i < fLength; i++)
		product *= *p++;

	return product;
}

/* max and min */
template <class nTYPE> 
nTYPE nArrayT<nTYPE>::Max(void) const
{
#if __option(extended_errorcheck)
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
	nTYPE  max   = *pthis++;
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
	nTYPE  max   = *pthis++;
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
	nTYPE  min   = *pthis++;
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE* pthis = Pointer();
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
	if (!fArray) throw eGeneralFail;
#endif

	/* ignore negative numbers */
	if (positive_only)
	{
		nTYPE* pthis = Pointer();
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
		nTYPE* pthis = Pointer();
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
	if (!fArray) throw eGeneralFail;
#endif

	nTYPE abs;
	nTYPE* pthis = Pointer();
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
nTYPE Dot(const nArrayT<nTYPE>& A1, const nArrayT<nTYPE>& A2)
{
/* dimension check */
#if __option (extended_errorcheck)
	if (A1.Length() != A2.Length()) throw eSizeMismatch;
#endif

	nTYPE* p1 = A1.Pointer();
	nTYPE* p2 = A2.Pointer();
	
	register nTYPE dot = 0.0;
	register nTYPE temp;
	
	int length = A1.Length();
	for (int i = 0; i < length; i++)
	{
		temp  = (*p1++);
		temp *= (*p2++);
		
		dot += temp;
	}
	
	return dot;
}

/* commonly used operations;
 *
 *      this =  scale*RHS;
 *		this += scale*RHS;
 *
 *  Note: returns a reference to (*this) for chaining operations */		
template <class nTYPE> 
void nArrayT<nTYPE>::SetToScaled(const nTYPE& scale, const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();

	nTYPE temp;

	for (int i = 0; i < fLength; i++)
	{
		temp  = scale;
		temp *= *pRHS++;	//multi-step needed incase pthis == pRHS

		*pthis++ = temp;
	}
}

template <class nTYPE> 
void nArrayT<nTYPE>::AddScaled(const nTYPE& scale, const nArrayT& RHS)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != RHS.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pRHS  = RHS.Pointer();
	
	nTYPE temp;
	
	for (int i = 0; i < fLength; i++)
	{
		temp  = scale;
		temp *= *pRHS++;
	
		*pthis++ += temp;
	}
}

/* sum (A + B) and difference (A - B) of nArrayT's */
template <class nTYPE> 
void nArrayT<nTYPE>::SumOf(const nArrayT& A, const nArrayT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength || fLength != B.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pA    = A.Pointer();
	nTYPE* pB    = B.Pointer();

	for (int i = 0; i < fLength; i++)
	{
		*pthis  = *pA++;
		*pthis += *pB++;
		
		pthis++;
	}
}
	
template <class nTYPE> 
void nArrayT<nTYPE>::DiffOf(const nArrayT& A, const nArrayT& B)
{
/* dimension checks */
#if __option (extended_errorcheck)
	if (fLength != A.fLength || fLength != B.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pA    = A.Pointer();
	nTYPE* pB    = B.Pointer();

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
	if (fLength != A.fLength || fLength != B.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pA    = A.Pointer();
	nTYPE* pB    = B.Pointer();

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
	    fLength != C.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pA    = A.Pointer();
	nTYPE* pB    = B.Pointer();
	nTYPE* pC    = C.Pointer();

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
	if (fLength != A.fLength || fLength != B.fLength) throw eSizeMismatch;
#endif

	nTYPE* pthis = Pointer();
	nTYPE* pA    = A.Pointer();
	nTYPE* pB    = B.Pointer();

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
void nArrayT<nTYPE>::AddCombination(const nArrayT& a, 
	const ArrayT<nArrayT<nTYPE>*>& A)
{
/* dimension of sum */
#if __option (extended_errorcheck)
	if (a.Length() != A.Length()) throw eSizeMismatch;
#endif

	/* sum */
	for (int i = 0; i < a.Length(); i++)
	{
		/* check each term in sum */
		if (fLength != A[i]->Length()) throw eSizeMismatch;
	
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
void nArrayT<nTYPE>::PrintWithFormat(ostream& out, int width, int prec,
	int wrapat, int tab) const
{
	int currprec = out.precision();
	out.precision(prec);
	
	int     count = 0;
	nTYPE* pthis = Pointer();
	
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

#endif /* _NARRAY_T_H_ */
