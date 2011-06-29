/* $Id: pMatrixT.h,v 1.3 2002-07-06 01:57:15 paklein Exp $ */

#ifndef _P_MATRIX_T_H_
#define _P_MATRIX_T_H_

/* base class */
#include "pArrayT.h"

namespace Tahoe {

/** A class to help working with a matrix of pointers. See pArrayT
 * for details of how the class handles the pointers. */
template <class TYPEPtr>
class pMatrixT: public pArrayT<TYPEPtr>
{
public:

	/** construct empty matrix */
	pMatrixT(void);

	/** construct matrix of pointers and initialize them all to NULL */
	pMatrixT(int numrows, int numcols);

	/** construct square matrix of pointers and initialize them all 
	 * to NULL */
	explicit pMatrixT(int squaredim);

	/** destructor */
	~pMatrixT(void);

	/** dimension matrix. Calls delete for every element of the matrix
	 * and initializes the entries in the new matrix to NULL */
	void Dimension(int numrows, int numcols);

	/** dimension matrix. Calls delete for every element of the matrix
	 * and initializes the entries in the new matrix to NULL */
	void Dimension(int squaredim);

	/** element accessor */
	ProxyTYPEPtr<TYPEPtr> operator()(int row, int col);

	/** const element accessor */
	ProxyTYPEPtr<TYPEPtr> operator()(int row, int col) const;

	/* dimensions */
	int Rows(void) const { return fRows; };
	int Cols(void) const { return fCols; };

private:

	/** copy/assignment operators */
	pMatrixT<TYPEPtr>& operator=(const pMatrixT& RHS);

protected:
	
	int	fRows;
	int	fCols;
};


/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPEPtr>
inline pMatrixT<TYPEPtr>::pMatrixT(void): fRows(0), fCols(0) { }

template <class TYPEPtr>
inline pMatrixT<TYPEPtr>::pMatrixT(int numrows, int numcols)
{
	Dimension(numrows, numcols);
}

template <class TYPEPtr>
inline pMatrixT<TYPEPtr>::pMatrixT(int squaredim)
{
	Dimension(squaredim);
}

/* destructor*/
template <class TYPEPtr>
inline pMatrixT<TYPEPtr>::~pMatrixT(void)
{
	fRows = 0;
	fCols = 0;
}

/* post construction dimensioning */
template <class TYPEPtr>
inline void pMatrixT<TYPEPtr>::Dimension(int numrows, int numcols)
{
	/* inherited */
	pArrayT<TYPEPtr>::Dimension(numrows*numcols);

	/* set dimensions */
	fRows = numrows;
	fCols = numcols;
}

template <class TYPEPtr>
inline void pMatrixT<TYPEPtr>::Dimension(int squaredim)
{
	Dimension(squaredim, squaredim);
}

/* element accessor */
template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pMatrixT<TYPEPtr>::operator()(int nrow, int ncol)
{
/* range checking */
#if __option (extended_errorcheck)
	if (nrow < 0 ||
	    nrow >= fRows ||
	    ncol < 0 ||
	    ncol >= fCols) throw eOutOfRange;
#endif
	
	return pArrayT<TYPEPtr>::operator[](ncol*fRows + nrow);
}

template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pMatrixT<TYPEPtr>::operator()(int nrow, int ncol) const
{
/* const_cast<> not supported */
#ifdef __SUNPRO_CC
	pArrayT<TYPEPtr>* const const_this = (pArrayT<TYPEPtr>* const) this;
#else
	pArrayT<TYPEPtr>* const const_this = const_cast<pArrayT<TYPEPtr>*>(this);
#endif

	return const_this->operator[](ncol*fRows + nrow);
}

} // namespace Tahoe 
#endif /* _P_MATRIX_T_H_ */
