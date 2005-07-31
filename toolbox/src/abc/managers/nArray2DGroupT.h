/* $Id: nArray2DGroupT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/16/1998)                                          */
/* Class to manage a list of equally-size nArray2DT<>'s. Storage          */
/* is grouped and all arrays added with Register can be set to new        */
/* major dimensions using Dimension. (percentover > 0) sets aside         */
/* extra space at memory allocation so that every call to Dimension       */
/* does not result in memory swapping.                                    */
/* NOTE: all registered arrays will be shallow.                           */

#ifndef _NARRAY2D_GROUP_T_H_
#define _NARRAY2D_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

/* direct members */
#include "nArray2DT.h"

template <class TYPE>
class nArray2DGroupT: public MemoryGroupT<TYPE>
{
public:

	/* constructor */
	nArray2DGroupT(int headroom);
	nArray2DGroupT(int headroom, int minordim);

	/* add Array2DT to list of managed */
	void Register(nArray2DT<TYPE>& array);

	/* (re-) dimension all arrays */
	void Dimension(int majordim, int minordim);
	void SetMinorDimension(int minordim); // does not allocate
	void SetMajorDimension(int majordim, bool copy_in);
		
private:

	/* size parameters */
	int fMajorDim;
	int fMinorDim;
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
	if (fMinorDim < 0) throw eGeneralFail;
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
	if (fMinorDim < 0) throw eGeneralFail;
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

#endif /* _NARRAY2D_GROUP_T_H_ */
