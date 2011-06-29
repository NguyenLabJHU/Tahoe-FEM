/* $Id: nMatrixGroupT.h,v 1.5 2005-05-01 19:27:51 paklein Exp $ */
/* created: paklein (04/17/1998) */
#ifndef _MATHMATRIX_GROUP_T_H_
#define _MATHMATRIX_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

/* direct members */
#include "nMatrixT.h"

namespace Tahoe {

/** class to manage a list of equally-size nMatrixT<>'s. Storage
 * is grouped and all matrices added with Register can be set to new
 * length using Dimension. (percentover > 0) sets aside
 * extra space at memory allocation so that every call to Dimension
 * does not result in memory swapping.
 * \note all registered matrices will be shallow.
 */
template <class TYPE>
class nMatrixGroupT: public MemoryGroupT<TYPE>
{
public:

	/** constructor */
	nMatrixGroupT(int headroom, bool pool_memory);

	/** add an nMatrixT to list of managed matricies */
	void Register(nMatrixT<TYPE>& matrix);

	/** (re-) dimension all matrices */
	void Dimension(int rows, int cols);	
	
private:

	/** \name dimensions of managed matricies */
	int fRows;
	int fCols;		
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
nMatrixGroupT<TYPE>::nMatrixGroupT(int headroom, bool pool_memory):
	MemoryGroupT<TYPE>(headroom, pool_memory),
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
		if (blocksize > this->BlockSize()) this->SetBlockSize(blocksize, false);

		/* reset pointers and dimensions */
		for (int i = 0; i < this->fArrays.Length(); i++)
		{
			/* safe cast due to type filtering by Register */
			nMatrixT<TYPE>* pmatrix = (nMatrixT<TYPE>*) this->fArrays[i];
		
			/* reset parameters */
			pmatrix->Set(fRows, fCols, this->BlockPointer(i));
		}
	}
}

} // namespace Tahoe 
#endif /* _MATHMATRIX_GROUP_T_H_ */
