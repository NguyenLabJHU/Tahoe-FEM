/* $Id: nArray2DGroupT.h,v 1.3 2002-10-20 22:38:53 paklein Exp $ */
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
	nArray2DGroupT(int headroom);
	nArray2DGroupT(int headroom, int minordim);
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
nArray2DGroupT<TYPE>::nArray2DGroupT(int headroom):
	MemoryGroupT<TYPE>(headroom),
	fMajorDim(0),
	fMinorDim(0)
{

}

template <class TYPE>
nArray2DGroupT<TYPE>::nArray2DGroupT(int headroom, int minordim):
	MemoryGroupT<TYPE>(headroom),
	fMajorDim(0),
	fMinorDim(minordim)
{
	/* error check */
	if (fMinorDim < 0) throw ExceptionT::kGeneralFail;
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
	if (majordim != fMajorDim || minordim != fMinorDim)
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
	if (fMinorDim < 0) throw ExceptionT::kGeneralFail;
}

template <class TYPE>
void nArray2DGroupT<TYPE>::SetMajorDimension(int majordim, bool copy_in)
{
	if (majordim != fMajorDim)
	{
		fMajorDim = majordim;

		/* need more memory, no memory reduction criteria */
		int blocksize = fMajorDim*fMinorDim;
		if (blocksize > BlockSize()) SetBlockSize(blocksize, copy_in);

		/* reset pointers and dimensions */
		for (int i = 0; i < fArrays.Length(); i++)
		{
			/* safe cast due to type filtering by Register */
			nArray2DT<TYPE>* parray = (nArray2DT<TYPE>*) fArrays[i];
		
			/* reset parameters */
			parray->Set(fMajorDim, fMinorDim, BlockPointer(i));
		}
	}
}

} // namespace Tahoe 
#endif /* _NARRAY2D_GROUP_T_H_ */
