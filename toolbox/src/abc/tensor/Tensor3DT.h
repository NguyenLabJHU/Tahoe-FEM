/* $Id: Tensor3DT.h,v 1.8 2005-04-30 21:29:51 paklein Exp $ */
/* created PAK (05/23/97) */

#ifndef _TENSOR3D_H_
#define _TENSOR3D_H_

/* base class */
#include "TensorT.h"

/* direct members */
#include "nMatrixT.h"

namespace Tahoe {

/* forward declarations */
template <class MATHTYPE> class nArrayT;

template <class MATHTYPE>
class Tensor3DT: public TensorT<MATHTYPE>
{
  public:

	/* constructors */
	Tensor3DT(void);
	Tensor3DT(int dim0, int dim1, int dim2);
	Tensor3DT(const Tensor3DT& source);

	/* dimensioning */
	void Dimension(int dim0, int dim1, int dim2);

	/** \deprecated replaced by Tensor3DT::Dimension on 02/13/2002 */
	void Allocate(int dim0, int dim1, int dim2) { Dimension(dim0, dim1, dim2); };

  	/* assignment operators */
  	Tensor3DT<MATHTYPE>& operator=(const Tensor3DT& RHS);
  	Tensor3DT<MATHTYPE>& operator=(const MATHTYPE& value);
	
	/** \name element and subdimension accessors. */
	/*@{*/
	MATHTYPE& operator()(int dim0, int dim1, int dim2);
	MATHTYPE* operator()(int dim0, int dim1);
	MATHTYPE* operator()(int dim0);

	const MATHTYPE& operator()(int dim0, int dim1, int dim2) const;
	const MATHTYPE* operator()(int dim0, int dim1) const;
	const MATHTYPE* operator()(int dim0) const;
	/*@{*/

	/*
	 * Contract the t3dex index of t3 with the t2dex index of
	 * t2 and place the results in this.  It is an error to
	 * pass (*this) as t3.  All dimensions must match exactly.
	 */
	void ContractIndex(const Tensor3DT& t3, int t3dex, 
	                   const nMatrixT<MATHTYPE>& t2, int t2dex);

	/*
	 * Contract the t3dex index of this with t1 and place the result
	 * in t2. All dimensions must match exactly.
	 */
	void ContractIndex(int t3dex, const nArrayT<MATHTYPE>& t1,
				       nMatrixT<MATHTYPE>& t2);
			
  private:
  
  	/*
  	 * Optimized constraction routines - routines assume all dimensions
  	 * have already been checked.
  	 */
	void Fast1Contract3D2D(const Tensor3DT& t3, 
  	                       const nMatrixT<MATHTYPE>& t2);
 	void Gen1Contract3D2D(const Tensor3DT& t3, int t3dex,
  	                      const nMatrixT<MATHTYPE>& t2, int t2dex);
	

  	/*
  	 * Optimized constraction routines - routines assume all dimensions
  	 * have already been checked.
  	 */
  	void Fast1Contract3D1D(const nArrayT<MATHTYPE>& t1, 
  	                       nMatrixT<MATHTYPE>& t2);
  	void Gen1Contract3D1D(int t3dex,
  	                      const nArrayT<MATHTYPE>& t1, 
  	                      nMatrixT<MATHTYPE>& t2);	

  protected:

	/* offsets */
	int fOffset0;
	int fOffset1;

};

/*************************************************************************
 *
 * Implementation
 *
 *************************************************************************/

/*
 * Constructor
 */
template <class MATHTYPE>
inline Tensor3DT<MATHTYPE>::Tensor3DT(void): fOffset0(0), fOffset1(0) { }

template <class MATHTYPE>
inline Tensor3DT<MATHTYPE>::Tensor3DT(int dim0, int dim1, int dim2)
{
	Dimension(dim0, dim1, dim2);
}

template <class MATHTYPE>
inline Tensor3DT<MATHTYPE>::Tensor3DT(const Tensor3DT& source)
{
	*this = source;
}

/*
 * Post-constructor
 */
template <class MATHTYPE> 
void Tensor3DT<MATHTYPE>::Dimension(int dim0, int dim1, int dim2)
{
	/* base class allocate */
	TensorT<MATHTYPE>::Dimension(dim0*dim1*dim2, 3);

	/* dimensions */
	this->fDim[0] = dim0;
	this->fDim[1] = dim1;
	this->fDim[2] = dim2;

	/* sanity check */
	if ( this->fDim.Min() < 1 ) throw ExceptionT::kOutOfRange;

	/* offsets */
	fOffset0 = this->fDim[1]*this->fDim[2];
	fOffset1 = this->fDim[2];
}

/*
 * Assignment operators
 */
template <class MATHTYPE> 
inline Tensor3DT<MATHTYPE>& Tensor3DT<MATHTYPE>::operator=(const Tensor3DT& RHS)
{
	/* inherited */
	TensorT<MATHTYPE>::operator=(RHS);
	return (*this);
}

template <class MATHTYPE> 
inline Tensor3DT<MATHTYPE>& Tensor3DT<MATHTYPE>::operator=(const MATHTYPE& value)
{
	/* inherited */
	TensorT<MATHTYPE>::operator=(value);
	return (*this);
}

/*
 * element and sub-dimension accessors.
 */
template <class MATHTYPE>
inline MATHTYPE& Tensor3DT<MATHTYPE>::operator()(int dim0, int dim1, int dim2)
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0] ||
	    dim1 < 0 || dim1 >= this->fDim[1] ||
	    dim2 < 0 || dim2 >= this->fDim[2]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray[dim0*fOffset0 + dim1*fOffset1 + dim2]);
}
template <class MATHTYPE>
inline const MATHTYPE& Tensor3DT<MATHTYPE>::operator()(int dim0, int dim1, int dim2) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0] ||
	    dim1 < 0 || dim1 >= this->fDim[1] ||
	    dim2 < 0 || dim2 >= this->fDim[2]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray[dim0*fOffset0 + dim1*fOffset1 + dim2]);
}

template <class MATHTYPE>
inline MATHTYPE* Tensor3DT<MATHTYPE>::operator()(int dim0, int dim1)
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0] ||
	    dim1 < 0 || dim1 >= this->fDim[1]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray + dim0*fOffset0 + dim1*fOffset1);
}
template <class MATHTYPE>
inline const MATHTYPE* Tensor3DT<MATHTYPE>::operator()(int dim0, int dim1) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0] ||
	    dim1 < 0 || dim1 >= this->fDim[1]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray + dim0*fOffset0 + dim1*fOffset1);
}

template <class MATHTYPE>
inline MATHTYPE* Tensor3DT<MATHTYPE>::operator()(int dim0)
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray + dim0*fOffset0);
}
template <class MATHTYPE>
inline const MATHTYPE* Tensor3DT<MATHTYPE>::operator()(int dim0) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (dim0 < 0 || dim0 >= this->fDim[0]) throw ExceptionT::kOutOfRange;
#endif

	return (this->fArray + dim0*fOffset0);
}

/*
 * Contract the t3dex index of t3 with the t2dex index of
 * t2 and place the results in this.  It is an error to
 * pass (*this) as t3.  All dimensions must match exactly.
 */
template <class MATHTYPE> 
void Tensor3DT<MATHTYPE>::ContractIndex(const Tensor3DT& t3, int t3dex, 
	const nMatrixT<MATHTYPE>& t2, int t2dex)
{
/* range checking */
#if __option (extended_errorcheck)
	if (this == &t3 || 
	    t3dex < 0 || t3dex >= 3 ||
        t2dex < 0 || t2dex >= 2) throw ExceptionT::kOutOfRange;
	
	/* (t2,t3) and (this, t2) */
	if (t2dex == 0)
		if (t2.Rows() != t3.fDim[t3dex] || t2.Cols() != this->fDim[2])
			throw ExceptionT::kOutOfRange;
	else if (t2.Cols() != t3.fDim[t3dex] || t2.Rows() != this->fDim[2]) 
		throw ExceptionT::kOutOfRange;
	
	/* (this, t3) */
	if (t3dex == 2)
		if (this->fDim[0] != t3.fDim[0] || this->fDim[1] != t3.fDim[1])
			throw ExceptionT::kOutOfRange;
	else if (t3dex == 1)
		if (this->fDim[0] != t3.fDim[0] || this->fDim[1] != t3.fDim[2])
			throw ExceptionT::kOutOfRange;
	else if (this->fDim[0] != t3.fDim[1] || this->fDim[1] != t3.fDim[2])
		throw ExceptionT::kOutOfRange;
#endif

	/* call optimized contraction routines */
	if (t2dex == 0 && t3dex == 2)
		Fast1Contract3D2D(t3, t2);
	else
		Gen1Contract3D2D(t3, t3dex, t2, t2dex);
}	                         

/*
 * Contract the t3dex index of this with t1 and place the result
 * in t2. All dimensions must match exactly.
 */
template <class MATHTYPE>
void Tensor3DT<MATHTYPE>::ContractIndex(int t3dex, 
	const nArrayT<MATHTYPE>& t1, nMatrixT<MATHTYPE>& t2)
{
#if __option (extended_errorcheck)
	if (t3dex < 0 || t3dex >= 3 ||t1.Length() != this->fDim[t3dex])
		throw ExceptionT::kOutOfRange;
	
	if (t3dex == 0)
		if (t2.Rows() != this->fDim[1] || t2.Cols() != this->fDim[2])
			throw ExceptionT::kOutOfRange;
	else if (t3dex == 1)
		if (t2.Rows() != this->fDim[0] || t2.Cols() != this->fDim[2])
			throw ExceptionT::kOutOfRange;
	else if (t2.Rows() != this->fDim[0] && t2.Cols() != this->fDim[1])
		throw ExceptionT::kOutOfRange;
#endif

	/* call optimized contraction routines */
	if (t3dex == 2 && 0)
		Fast1Contract3D1D(t1, t2);
	else
		Gen1Contract3D1D(t3dex, t1, t2);
}	                         

/**********************************************************************
 *
 *  Private
 *
 **********************************************************************/

/*
 * Optimized constraction routines - routines assume all dimensions
 * have already been checked.
 */

/*
 * computes: this_ijk = t3_ijm t2_mk
 */ 
template <class MATHTYPE>
void Tensor3DT<MATHTYPE>::Fast1Contract3D2D(const Tensor3DT& t3, 
	const nMatrixT<MATHTYPE>& t2)
{
	/* increments */
	int t2kinc = t2.Rows();
	
	int t3iinc = t3.fOffset0;
	int t3jinc = t3.fOffset1;
	
	/* first ptr's */
	MATHTYPE* pData = this->Pointer();
	MATHTYPE* pt2   = t2.Pointer();
	MATHTYPE* pt3   = t3.Pointer();
	
	int dotrange = t3.fDim[2];

	register MATHTYPE sum;
	register MATHTYPE temp;

	for (int i = 0; i < this->fDim[0]; i++)
	{
		MATHTYPE* pt3j = pt3;
	
		for (int j = 0; j < this->fDim[1]; j++)
		{
			MATHTYPE* pt2k = pt2;
		
			for (int k = 0; k < this->fDim[2]; k++)
			{
				sum = 0.0;
			
				MATHTYPE* pt2m = pt2k;
				MATHTYPE* pt3m = pt3j;
			
				/* contraction */	
				for (int m = 0; m < dotrange; m++)
				{
					temp  = *pt2m++;
					temp *= *pt3m++;
					
					sum += temp;
				}

				*pData++ = sum;
								
				pt2k += t2kinc;
			}
			
			pt3j += t3jinc;
		}
		
		pt3 += t3iinc;
	}
}

/*
 * computes: this_ijk = t3_ijm t2_mk
 *           ...
 *           this_ijk = t3_mij t2_km
 */ 
template <class MATHTYPE>
void Tensor3DT<MATHTYPE>::Gen1Contract3D2D(const Tensor3DT& t3, int t3dex,
	const nMatrixT<MATHTYPE>& t2, int t2dex)
{
	/* set t2 increments */
	int t2minc;
	int t2kinc;
	if (t2dex == 0) {
		t2minc = 1;
		t2kinc = t2.Rows();
	} else {
		t2minc = t2.Rows();
		t2kinc = 1;
	}
	
	/* set t3 increments */
	int t3minc;
	int t3jinc;
	int t3iinc;
	if (t3dex == 2) {
		t3minc = 1;
		t3jinc = t3.fOffset1;
		t3iinc = t3.fOffset0;
	} else if (t3dex == 1) {
		t3minc = t3.fOffset1;
		t3jinc = 1;
		t3iinc = t3.fOffset0;
	} else {
		t3minc = t3.fOffset0;
		t3jinc = 1;
		t3iinc = t3.fOffset1;
	}
			
	/* first ptr's */
	MATHTYPE* pData = this->Pointer();
	MATHTYPE* pt3   = t3.Pointer();
	MATHTYPE* pt2   = t2.Pointer();
		
	int dotrange = t3.fDim[t3dex];

	register MATHTYPE sum;
	register MATHTYPE temp;
	
	for (int i = 0; i < this->fDim[0]; i++)
	{
		MATHTYPE* pt3j = pt3;
	
		for (int j = 0; j < this->fDim[1]; j++)
		{
			MATHTYPE* pt2k = pt2;
		
			for (int k = 0; k < this->fDim[2]; k++)
			{
				sum = 0.0;
			
				MATHTYPE* pt2m = pt2k;
				MATHTYPE* pt3m = pt3j;
			
				/* contraction */		
				for (int m = 0; m < dotrange; m++)
				{
					temp  = *pt2m;
					temp *= *pt3m;

					sum  += temp;
			
					pt2m += t2minc;
					pt3m += t3minc;
				}
				*pData++ = sum;
				
				pt2k += t2kinc;
			}
			
			pt3j += t3jinc;
		}
		
		pt3 += t3iinc;
	}	
}  	                        	

/*
 * Optimized constraction routines - routines assume all dimensions
 * have already been checked.
 */

/*
 * returns: t2_ij = this_ijk t1_k
 */ 
template <class MATHTYPE> 
void Tensor3DT<MATHTYPE>::Fast1Contract3D1D(const nArrayT<MATHTYPE>& t1, 
	nMatrixT<MATHTYPE>& t2)
{
	/* increments */
	int t2jinc = t2.Rows();
	int t3iinc = fOffset0;
	int t3jinc = fOffset1;
	
	/* first ptr's */
	MATHTYPE* pt2 = t2.Pointer();
	MATHTYPE* pt3 = this->Pointer();
	MATHTYPE* pt1 = t1.Pointer();
	
	int dotrange = this->fDim[2];
	
	register MATHTYPE sum;
	register MATHTYPE temp;
	
	for (int i = 0; i < this->fDim[0]; i++)
	{
		MATHTYPE* pt3j  = pt3;
		MATHTYPE* pData = pt2;
	
		for (int j = 0; j < this->fDim[1]; j++)
		{		
			sum = 0.0;
							
			MATHTYPE* pt1m = pt1;
			MATHTYPE* pt3m = pt3j;
			
			/* contraction */	
			for (int m = 0; m < dotrange; m++)
			{
				temp  = *pt1m++;
				temp *= *pt3m++;
			
				sum  += temp;
			}
							
			*pData = sum;
			
			pData += t2jinc;
			pt3j  += t3jinc;
		}
		
		pt2++;
		pt3 += t3iinc;
	}
}

/*
 * returns: t2_ij = this_ijm t1_m	-or-	
 *			      = this_imj t1_m   -or-
 *                = this_mij t1_m
 */ 
template <class MATHTYPE> 
void Tensor3DT<MATHTYPE>::Gen1Contract3D1D(int t3dex, 
	const nArrayT<MATHTYPE>& t1, 
	nMatrixT<MATHTYPE>& t2)
{
	/* set t3 increments */
	int t3minc;
	int t3jinc;
	int t3iinc;
	if (t3dex == 2) {
		t3minc = 1;
		t3jinc = fOffset1;
		t3iinc = fOffset0;
	} else if (t3dex == 1) {
		t3minc = fOffset1;
		t3jinc = 1;
		t3iinc = fOffset0;
	} else {
		t3minc = fOffset0;
		t3jinc = 1;
		t3iinc = fOffset1;
	}

	/* t2 increment */
	int t2jinc = t2.Rows();			
			
	/* first ptr's */
	MATHTYPE* pt2   = t2.Pointer();
	MATHTYPE* pt3   = this->Pointer();
	MATHTYPE* pt1   = t1.Pointer();
		
	int dotrange = this->fDim[t3dex];

	register MATHTYPE sum;
	register MATHTYPE temp;
	
	for (int i = 0; i < t2.Rows(); i++)
	{
		MATHTYPE* pt3j  = pt3;
		MATHTYPE* pData = pt2; 
	
		for (int j = 0; j < t2.Cols(); j++)
		{
			MATHTYPE* pt1m = pt1;
			MATHTYPE* pt3m = pt3j;
		
			sum = 0.0;			
			
			/* contraction */		
			for (int m = 0; m < dotrange; m++)
			{
				temp  = *pt1m++;
				temp *= *pt3m;
			
				sum += temp;
			
				pt3m += t3minc;
			}

			*pData = sum;
				
			pData += t2jinc;			
			pt3j  += t3jinc;
		}
		
		pt3 += t3iinc;
		pt2++;
	}	
}  	                        	

}//namespace Tahoe
#endif /* _TENSOR3D_H_ */
