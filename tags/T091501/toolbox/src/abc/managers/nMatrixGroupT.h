/* $Id: nMatrixGroupT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/17/1998)                                          */
/* Class to manage a list of equally-size nMatrixT<>'s. Storage           */
/* is grouped and all matrices added with Register can be set to new      */
/* length using Dimension. (percentover > 0) sets aside                   */
/* extra space at memory allocation so that every call to Dimension       */
/* does not result in memory swapping.                                    */
/* NOTE: all registered matrices will be shallow.                         */

#ifndef _MATHMATRIX_GROUP_T_H_
#define _MATHMATRIX_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

/* direct members */
#include "nMatrixT.h"

template <class TYPE>
class nMatrixGroupT: public MemoryGroupT<TYPE>
{
public:

	/* constructor */
	nMatrixGroupT(int headroom);

	/* add Array2DT to list of managed */
	void Register(nMatrixT<TYPE>& matrix);

	/* (re-) dimension all matrices */
	void Dimension(int rows, int cols);	
	
private:

	int fRows;
	int fCols;		
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
nMatrixGroupT<TYPE>::nMatrixGroupT(int headroom):
	MemoryGroupT<TYPE>(headroom),
	fRows(0),
	fCols(0)
{

}

/* add nMatrixT to list of managed - function allows only nMatrixT's
* to be registered, ie. an matrix type filter */
template <class TYPE>
inline void nMatrixGroupT<TYPE>::Register(nMatrixT<TYPE>& matrix)
{
	/* inherited */
	MemoryGroupT<TYPE>::Register(matrix);
}

/* (re-) dimension all matrices */
template <class TYPE>
void nMatrixGroupT<TYPE>::Dimension(int rows, int cols)
{
	if (rows != fRows || cols != fCols)
	{
		fRows = rows;
		fCols = cols;
	
		/* need more memory, no memory reduction criteria */
		int blocksize = fRows*fCols;
		if (blocksize > BlockSize()) SetBlockSize(blocksize, false);

		/* reset pointers and dimensions */
		for (int i = 0; i < fArrays.Length(); i++)
		{
			/* safe cast due to type filtering by Register */
			nMatrixT<TYPE>* pmatrix = (nMatrixT<TYPE>*) fArrays[i];
		
			/* reset parameters */
			pmatrix->Set(fRows, fCols, BlockPointer(i));
		}
	}
}

#endif /* _MATHMATRIX_GROUP_T_H_ */
