/* $Id: nArrayGroupT.h,v 1.6 2005-06-05 06:21:07 paklein Exp $ */
/* created: paklein (04/17/1998) */
#ifndef _NARRAY_GROUP_T_H_
#define _NARRAY_GROUP_T_H_

/* base class */
#include "MemoryGroupT.h"

namespace Tahoe {

/** class to manage a list of equally-size nArrayT<>'s. Storage
 * is grouped and all arrays added with Register can be set to new 
 * length using Dimension. (percentover > 0) sets aside extra space
 * at memory allocation so that every call to nArrayGroupT::Dimension 
 * does not result in memory swapping.
 * \note all registered arrays will be shallow.
 */
template <class TYPE>
class nArrayGroupT: public MemoryGroupT<TYPE>
{
public:

	/** constructor */
	nArrayGroupT(int headroom, bool pool_memory);

	/** add an nArrayT to list of managed arrays */
	int Register(nArrayT<TYPE>& array);

	/** access a registered array */
	nArrayT<TYPE>& Array(int i);

	/** (re-)dimension all arrays and copy in old data if specified */
	void Dimension(int length, bool copy_in);
	
private:

	/** current length of group set by call to nArrayGroupT::Length */
	int fLength;
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructor */
template <class TYPE>
nArrayGroupT<TYPE>::nArrayGroupT(int headroom, bool pool_memory):
	MemoryGroupT<TYPE>(headroom, pool_memory),
	fLength(0)
{

}

/* access a registered array */
template <class TYPE>
nArrayT<TYPE>& nArrayGroupT<TYPE>::Array(int i)
{
	/* cast is OK because Register accepts only nArrayT's */
	nArrayT<TYPE>* p = (nArrayT<TYPE>*) this->fArrays[i];
	return *p;
}

/* add nArrayT to list of managed - function allows only nArrayT's
* to be registered, ie. an array type filter */
template <class TYPE>
inline int nArrayGroupT<TYPE>::Register(nArrayT<TYPE>& array)
{
	/* inherited */
	int index = MemoryGroupT<TYPE>::Register(array);
	
	/* return dimensioned array */
	if (fLength > 0)
		array.Alias(fLength, this->BlockPointer(index));
	else /* empty array */
		array.Dimension(0);
	
	return index;
}

/* (re-) dimension all arrays */
template <class TYPE>
void nArrayGroupT<TYPE>::Dimension(int length, bool copy_in)
{
	if (fLength != length)
	{
		fLength = length;
	
		/* need more memory, no memory reduction criteria */
		if (fLength > this->BlockSize()) this->SetBlockSize(fLength, copy_in);

		/* reset pointers and dimensions */
		for (int i = 0; i < this->fArrays.Length(); i++)
			this->fArrays[i]->Alias(fLength, this->BlockPointer(i));
	}
}

} // namespace Tahoe 
#endif /* _NARRAY_GROUP_T_H_ */
