/* $Id: TensorT.h,v 1.4 2005-04-30 21:14:36 paklein Exp $ */
/* created paklein (05/19/97) */
#ifndef _TENSORBASET_H_
#define _TENSORBASET_H_

/* ANSI headers */
#include <iostream.h>

/* base class */
#include "nArrayT.h"

/* direct members */
#include "iArrayT.h"


namespace Tahoe {

template <class MATHTYPE>
class TensorT: public nArrayT<MATHTYPE>
{
  public:

	/* constructor */
	TensorT(void);
	TensorT(int size, int rank);
	TensorT(const TensorT& source);

	/* set dimension */
	void Dimension(int size, int rank);

	/* set dimension */
	void Allocate(int size, int rank) { Dimension(size, rank); };

	/* dimensions */
	int Rank(void) const; 
	int Dim(int dim) const;

  	/* copy/assignment operators */
  	TensorT<MATHTYPE>& operator=(const TensorT& RHS);
  	TensorT<MATHTYPE>& operator=(const MATHTYPE& value);
  	void ShallowCopy(const TensorT& RHS);

	/* print parsed tensor to the ostream */
	void PrintTensor(ostream& out) const;

  protected:
  	
  protected:
	
	/* dimensions */
	iArrayT	fDim;
};

/*************************************************************************
 *
 * Implementation
 *
 *************************************************************************/

/*
 * Constructors
 */
template <class MATHTYPE> 
inline TensorT<MATHTYPE>::TensorT(void) { }

template <class MATHTYPE> 
inline TensorT<MATHTYPE>::TensorT(int size, int rank)
{
	Dimension(size, rank);	
}

template <class MATHTYPE> 
inline TensorT<MATHTYPE>::TensorT(const TensorT& source)
{
	*this = source;
}

/*
 * Post-constructor
 */
template <class MATHTYPE> 
void TensorT<MATHTYPE>::Dimension(int size, int rank)
{
	/* dimensions */
	fDim.Dimension(rank);
	
	/* base class */
	nArrayT<MATHTYPE>::Dimension(size);
}

/*
 * Dimensions
 */
template <class MATHTYPE> 
inline int TensorT<MATHTYPE>::Rank(void) const
{
	return ( fDim.Length() );
}

template <class MATHTYPE> 
inline int TensorT<MATHTYPE>::Dim(int dim) const
{
	return (fDim[dim]);
}

/*
 * Copy/assignment operators
 */
template <class MATHTYPE> 
inline TensorT<MATHTYPE>& TensorT<MATHTYPE>::operator=(const TensorT& RHS)
{
	/* inherited */
	nArrayT<MATHTYPE>::operator=(RHS);

	/* copy dimensions */
	fDim = RHS.fDim;

	return (*this);
}

template <class MATHTYPE> 
inline TensorT<MATHTYPE>& TensorT<MATHTYPE>::operator=(const MATHTYPE& value)
{
	/* inherited */
	nArrayT<MATHTYPE>::operator=(value);
	return (*this);
}

template <class MATHTYPE> 
inline void TensorT<MATHTYPE>::ShallowCopy(const TensorT& RHS)
{
	/* (deep) copy dimensions */
	fDim = RHS.fDim;

	/* base class */
	nArrayT<MATHTYPE>::ShallowCopy(RHS);
}

/*
 * Input operator
 */
template <class MATHTYPE>
istream& operator>>(istream& in, TensorT<MATHTYPE>& RHS)
{
	MATHTYPE* p = RHS.Pointer();
	for (int i = 0; i < RHS.Length(); i++)
		in >> *p++;
	return (in);
}

/*
 * Print parsed tensor to the ostream
 */
template <class MATHTYPE> 
void TensorT<MATHTYPE>::PrintTensor(ostream& out) const
{
	iArrayT dex( Rank() );
	dex = 0;

	MATHTYPE* pD = this->fData;
	for (int i = 0; i < this->fTotalSize; i++)
	{
		/* print indices */
		for (int j = 0; j < this->fRank; j++)
			out << setw(kIntWidth) << dex[j];
		out << *pD++ << '\n';

		/* increment indices */
		dex[this->fRank-1]++;
		for (int j = this->fRank-1; j > 0; j--)
			if (dex[j] == fDim[j])
			{
				dex[j] = 0;
				dex[j-1]++;
			}
	}

	out << '\n';
}

/**********************************************************************
 *  Protected
 **********************************************************************/

/*
 * Rank and dimension checking function.  Returns 1 if the rank and
 * all the dimensions are identical and returns 0 otherwise.
 */
template <class MATHTYPE1, class MATHTYPE2> 
int SameDimensions(const TensorT<MATHTYPE1>& t1, const TensorT<MATHTYPE2>& t2)
{
	if (t1.Rank() != t2.Rank())
		return (0);
	
	for (int i = 0; i < t1.Rank(); i++)
		if (t1.Dim(i) != t2.Dim(i))
			return (0);

	return (1);
}

} // namespace Tahoe 
#endif /* _TENSORBASET_H_ */
