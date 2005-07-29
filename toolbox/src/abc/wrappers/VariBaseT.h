/* $Id: VariBaseT.h,v 1.7 2005-07-29 03:09:35 paklein Exp $ */
/* created: paklein (04/18/1998) */
#ifndef _VARI_BASE_T_H_
#define _VARI_BASE_T_H_

/* direct members */
#include "ArrayT.h"

namespace Tahoe {

/** \name base class for WRAPPERS of ArrayT<>'s, and derivatives.
 * Adds dynamic re-sizing using some headroom to cut down
 * calls for memory de/re-allocation */
template <class TYPE>
class VariBaseT
{
public:

	/** \name constructors */
	/*@{*/
	VariBaseT(void);
	
	/** construct and define headroom
	 * \param headroom amount in percent of excess memory that is allocated */
	VariBaseT(int headroom);
	/*@}*/

	/* set the head room parameter
	 * \param headroom amount in percent of excess memory that is allocated */
	void SetHeadRoom(int headroom);

protected:
	
	/** \name convert array
	 * convert array to shallow copy of fMemory, fill
	 * extra new space and copy old data, if specified */
	/*@{*/
	void SetAlias(ArrayT<TYPE>& array, int length, bool copy_in);
	void SetAlias(ArrayT<TYPE>& array, int length, const TYPE& fill, bool copy_in);
	/*@}*/

	/** free memory */
	void Free(void);

	/** swap memory with the source array */
	void Swap(ArrayT<TYPE>& source);

private:

	/** \name not allowed */
	/*@{*/
	/** no copy constructor */
	VariBaseT(const VariBaseT& source);

	/** no assigment operator */	 			  	
	void operator=(const VariBaseT& RHS);
	/*@}*/

private:

	/** overhead size - % */
	int fHeadRoom;
	
	/** memory space */
	ArrayT<TYPE> fMemory;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
VariBaseT<TYPE>::VariBaseT(void): fHeadRoom(0) { }
template <class TYPE>
VariBaseT<TYPE>::VariBaseT(int headroom) { SetHeadRoom(headroom); }

/* set the head room parameter */
template <class TYPE>
inline void VariBaseT<TYPE>::SetHeadRoom(int headroom)
{
	fHeadRoom = headroom;
	if (fHeadRoom < 0) ExceptionT::GeneralFail("VariBaseT");
}

/**********************************************************************
*  Protected
**********************************************************************/
	
/* set length of the ward, fill extra space if specified */
template <class TYPE>
void VariBaseT<TYPE>::SetAlias(ArrayT<TYPE>& array, int length, bool copy_in)
{
	if (array.Length() != length || array.IsAllocated())
	{
		/* need more memory (no criteria to reallocate smaller) */
		if (length > fMemory.Length())
		{
			int memsize = (length*(100 + fHeadRoom))/100;
			
			/* first time */
			if (array.IsAllocated())	
			{
				/* allocate space */
				fMemory.Dimension(memsize);
		
				/* copy data from the ward */
				if (copy_in)
				{
					int copysize = (array.Length() < length) ? array.Length():length;
					fMemory.CopyPart(0, array, 0, copysize);
				}
			}
			else if (copy_in)
				fMemory.Resize(memsize);
			else
				fMemory.Dimension(memsize);
		}
			
		/* ward becomes shallow copy */
		array.Set(length, fMemory.Pointer());
	}
}

template <class TYPE>
void VariBaseT<TYPE>::SetAlias(ArrayT<TYPE>& array, int length,
	const TYPE& fill, bool copy_in)
{
	/* get old length */
	int oldlength = array.Length();

	/* resize */
	SetAlias(array, length, copy_in);
	
	/* initialize added space */
	TYPE* parray = array.Pointer() + oldlength;		
	for (int i = oldlength; i < length; i++)
		*parray++ = fill;
}

/* free memory */
template <class TYPE>
inline void VariBaseT<TYPE>::Free(void) { fMemory.Free(); }

/* swap memory */
template <class TYPE>
inline void VariBaseT<TYPE>::Swap(ArrayT<TYPE>& source) { fMemory.Swap(source); }

} // namespace Tahoe 
#endif /* _VARI_BASE_T_H_ */
