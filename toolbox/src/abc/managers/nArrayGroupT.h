/* $Id: nArrayGroupT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/17/1998)                                          */
/* Class to manage a list of equally-size nArrayT<>'s. Storage            */
/* is grouped and all arrays added with Register can be set to new        */
/* length using Dimension. (percentover > 0) sets aside                   */
/* extra space at memory allocation so that every call to Dimension       */
/* does not result in memory swapping.                                    */
/* NOTE: all registered arrays will be shallow.                           */

#ifndef _NARRAY_GROUP_T_H_
#define _NARRAY_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

template <class TYPE>
class nArrayGroupT: public MemoryGroupT<TYPE>
{
public:

	/* constructor */
	nArrayGroupT(int headroom);

	/* add Array2DT to list of managed */
	void Register(nArrayT<TYPE>& array);

	/* (re-) dimension all arrays and copy in old data if specified */
	void Dimension(int length, bool copy_in);
	
private:

	int fLength; // current group length		
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructor */
template <class TYPE>
nArrayGroupT<TYPE>::nArrayGroupT(int headroom):
	MemoryGroupT<TYPE>(headroom),
	fLength(0)
{

}

/* add nArrayT to list of managed - function allows only nArrayT's
* to be registered, ie. an array type filter */
template <class TYPE>
inline void nArrayGroupT<TYPE>::Register(nArrayT<TYPE>& array)
{
	/* inherited */
	MemoryGroupT<TYPE>::Register(array);
}

/* (re-) dimension all arrays */
template <class TYPE>
void nArrayGroupT<TYPE>::Dimension(int length, bool copy_in)
{
	if (fLength != length)
	{
		fLength = length;
	
		/* need more memory, no memory reduction criteria */
		if (fLength > BlockSize()) SetBlockSize(fLength, copy_in);

		/* reset pointers and dimensions */
		for (int i = 0; i < fArrays.Length(); i++)
			fArrays[i]->Set(fLength, BlockPointer(i));
	}
}

#endif /* _NARRAY_GROUP_T_H_ */
