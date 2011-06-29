/* $Id: nArray2DGroupT.h,v 1.6 2005-05-01 19:27:51 paklein Exp $ */
/* created: paklein (04/16/1998) */
#ifndef _NARRAY2D_GROUP_T_H_
#define _NARRAY2D_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

/* direct members */
#include "nArray2DT.h"

namespace Tahoe {

/** class to manage a list of equally-size nArray2DT<>'s. Storage
 * is grouped and all arrays added with Register can be set to new
 * major dimensions using Dimension. (percentover > 0) sets aside
 * extra space at memory allocation so that every call to nArray2DGroupT::Dimension
 * does not result in memory swapping.
 * \note all registered arrays will be shallow.
 */
template <class TYPE>
class nArray2DGroupT: public MemoryGroupT<TYPE>
{
public:

	/** \name constructors */
	/*@{*/
	nArray2DGroupT(int headroom, bool pool_memory);
	nArray2DGroupT(int headroom, bool pool_memory, int minordim);
	/*@}*/

	/** add an nArray2DT to list of managed arrays */
	void Register(nArray2DT<TYPE>& array);

	/** (re-)dimension all arrays */
	/*@{*/
	void Dimension(int majordim, int minordim);
	void SetMajorDimension(int majordim, bool copy_in);

	/** set minor dimension of but does not reset the managed
	 * arrays. Arrays are dimensioned by the next call to
	 * nArray2DGroupT::SetMajorDimension. */
	void SetMinorDimension(int minordim);
	/*@}*/
		
private:

	/** \name size parameters */
	/*@{*/
	int fMajorDim;
	int fMinorDim;
	/*@}*/
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructors */
template <class TYPE>
nArray2DGroupT<TYPE>::nArray2DGroupT(int headroom, bool pool_memory):
	MemoryGroupT<TYPE>(headroom, pool_memory),
	fMajorDim(0),
	fMinorDim(0)
{

}

template <class TYPE>
nArray2DGroupT<TYPE>::nArray2DGroupT(int headroom, bool pool_memory, int minordim):
	MemoryGroupT<TYPE>(headroom, pool_memory),
	fMajorDim(0),
	fMinorDim(minordim)
{
	/* error check */
	if (fMinorDim < 0)
		ExceptionT::GeneralFail("nArray2DGroupT<TYPE>::nArray2DGroupT",
			"bad dimensiion %d", fMinorDim);
}

/* add Array2DT to list of managed - function allows only nArray2DT's
* to be registered, ie. an array type filter */
template <class TYPE>
inline void nArray2DGroupT<TYPE>::Register(nArray2DT<TYPE>& array)
{
	/* inherited */
	MemoryGroupT<TYPE>::Register(array);
}

/* (re-) dimension all arrays */
template <class TYPE>
inline void nArray2DGroupT<TYPE>::Dimension(int majordim, int minordim)
{
	if (majordim != fMajorDim || minordim != fMinorDim || fMajorDim == 0) /* 0 => uninitialized */
	{
		/* reset dimensions */
		fMinorDim = minordim;
		fMajorDim = majordim - 1; // ensure SetMajorDimension reallocates
		
		/* set rest (don't copy old data) */
		SetMajorDimension(majordim, false);
	}
}

template <class TYPE>
void nArray2DGroupT<TYPE>::SetMinorDimension(int minordim)
{
	fMinorDim = minordim;

	/* error check */
	if (fMinorDim < 0)
		ExceptionT::GeneralFail("nArray2DGroupT<TYPE>::SetMinorDimension", 
			"bad dimension %d", fMinorDim);
}

template <class TYPE>
void nArray2DGroupT<TYPE>::SetMajorDimension(int majordim, bool copy_in)
{
	if (majordim != fMajorDim || fMajorDim == 0) /* 0 => uninitialized */
	{
		fMajorDim = majordim;

		/* need more memory, no memory reduction criteria */
		int blocksize = fMajorDim*fMinorDim;
		if (blocksize > this->BlockSize()) this->SetBlockSize(blocksize, copy_in);

		/* reset pointers and dimensions */
		for (int i = 0; i < this->fArrays.Length(); i++)
		{
			/* safe cast due to type filtering by Register */
			nArray2DT<TYPE>* parray = (nArray2DT<TYPE>*) this->fArrays[i];
		
			/* reset parameters */
			parray->Alias(fMajorDim, fMinorDim, this->BlockPointer(i));
		}
	}
}

} /* namespace Tahoe */

#endif /* _NARRAY2D_GROUP_T_H_ */
