/* $Id: AutoFill2DT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (01/19/1999)                                          */
/* NOTE: going to use this with a manager to help count and               */
/* store edges in a graph. There have to be 2 modes of                    */
/* usage.                                                                 */
/* (1) just let AutoFill2DT go nuts with managing itself                  */
/* (2) try to implement a "shuffle down", where vacated spots             */
/* in the array are filled with same size blocks pulled                   */
/* from the end of the list.                                              */
/* Array that automatically increases its minor dimension when            */
/* elements are inserted into the rows using Append() or                  */
/* AppendUnique().                                                        */
/* NOTE: An overloaded assignment operator is not used for                */
/* duplicating the array contents. Contents is simply                     */
/* byte copied.                                                           */

#ifndef _AUTO_ARRAY2D_T_H_
#define _AUTO_ARRAY2D_T_H_

#include <string.h>

#include "Environment.h"
#include "ExceptionCodes.h"

template <class TYPE>
class AutoFill2DT
{
public:

	/* constructors */
	AutoFill2DT(int majordim, int headroom);
	AutoFill2DT(int majordim, int headroom, int maxminordim);
	AutoFill2DT(const AutoFill2DT<TYPE>& source);

	/* destructor */
	~AutoFill2DT(void);

	/* dimensions */
	int MajorDim(void) const;
	int HeadRoom(void) const;
	int MinorDim(int row) const;
	int MaxMinorDim(void) const;
	int MinMinorDim(void) const;	
	int LogicalSize(void) const; // size of data (<= total memory allocation)
	int MinorDimCount(int dim) const; // returns number of rows with the specified size

	/* re-sizing */
	void SetMaxMinorDim(int maxminordim, int make_headroom = 0);
	void SetHeadRoom(int headroom);
	
	/* accessors */
	TYPE& operator()(int majordim, int minordim) const;
	TYPE* operator()(int majordim) const;

	/* set logical size = 0 */
	void Reset(void);         // for all rows
	void Reset(int majordim); // for the selected row

	/* copy of source - dimensions must match */
	void Copy(const AutoFill2DT<TYPE>& source);
	
	/* adding values to rows */
	void Append(int majordim, const TYPE& value);
	void Append(int majordim, const ArrayT<TYPE>& source);
	
	int AppendUnique(int majordim, const TYPE& value);          // returns 1 if value added
	int AppendUnique(int majordim, const ArrayT<TYPE>& source); // returns number appended
	
private:

	/* resize data and copy in previous values. fMaxMinorDim grows by at
	 * least 1 with every call */
	void Allocate(int maxminordim, int headroom);
	
	/* no assigment operator */
	AutoFill2DT& operator=(AutoFill2DT&);
	
private:

	/* dimensions */
	int fMajorDim;
	int fHeadRoom; // percent over-allocation
	int fMaxMinorDim;

	/* row counts */
	int* fCounts;

	/* the data */
	TYPE* fArray;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
inline AutoFill2DT<TYPE>::AutoFill2DT(int majordim, int headroom):
	fMajorDim(majordim),
	fHeadRoom(headroom),
	fMaxMinorDim(0),
	fCounts(NULL),
	fArray(NULL)	
{
	/* check */
	if (fMajorDim < 0 ||
	    fHeadRoom < 0) throw eGeneralFail;	
}

template <class TYPE>
//inline
AutoFill2DT<TYPE>::AutoFill2DT(int majordim, int headroom, int maxminordim):
	fMajorDim(majordim),
	fHeadRoom(headroom),
	fMaxMinorDim(0),
	fCounts(NULL),
	fArray(NULL)	
{
	/* check */
	if (fMajorDim < 0 ||
		fMaxMinorDim < 0 ||
	    fHeadRoom < 0) throw eGeneralFail;	

	/* initialize memory */
	Allocate(maxminordim, 0);
}

template <class TYPE>
AutoFill2DT<TYPE>::AutoFill2DT(const AutoFill2DT<TYPE>& source):
	fMajorDim(source.fMajorDim),
	fHeadRoom(source.fHeadRoom),
	fMaxMinorDim(0),
	fCounts(NULL),
	fArray(NULL)	
{
	/* allocate memory */
	Allocate(source.fMaxMinorDim, 0);

	/* copy data */
	this->Copy(source);
}

template <class TYPE>
AutoFill2DT<TYPE>::~AutoFill2DT(void)
{
	delete[] fCounts;
	fCounts = NULL;

	delete[] fArray;
	fArray = NULL;
}

/* dimensions */
template <class TYPE>
inline int AutoFill2DT<TYPE>::MajorDim(void) const { return fMajorDim; }

template <class TYPE>
int AutoFill2DT<TYPE>::HeadRoom(void) const { return fHeadRoom; }

template <class TYPE>
inline int AutoFill2DT<TYPE>::MinorDim(int majordim) const
{
#if __option(extended_errorcheck)
	/* range check */
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	return fCounts[majordim];
}

template <class TYPE>
int AutoFill2DT<TYPE>::MaxMinorDim(void) const { return fMaxMinorDim; }

template <class TYPE>
int AutoFill2DT<TYPE>::MinMinorDim(void) const
{
	int  min    = *fCounts;
	int* pcount =  fCounts + 1;

	for (int i = 1; i < fMajorDim; i++)
	{
		min  = (*pcount < min) ? *pcount : min;
		pcount++;
	}
	
	return min;
}

template <class TYPE>
int AutoFill2DT<TYPE>::LogicalSize(void) const
{
	int    size = 0;
	int* pcount = fCounts;

	for (int i = 0; i < fMajorDim; i++)
		size += *pcount++;
	
	return size;
}

template <class TYPE>
int AutoFill2DT<TYPE>::MinorDimCount(int dim) const // returns number of rows with the specified size
{
	int   count = 0;
	int* pcount = fCounts;
	for (int i = 0; i < fMajorDim; i++)
		if (*pcount++ == dim) count++;

	return count;
}

/* re-sizing */
template <class TYPE>
inline void AutoFill2DT<TYPE>::SetMaxMinorDim(int maxminordim, int make_headroom)
{
	//TEMP - size can only grow for now
	if (maxminordim < fMaxMinorDim) throw eGeneralFail;

	/* set memory */
	Allocate(maxminordim, (make_headroom == 1) ? fHeadRoom : 0);
}

template <class TYPE>
inline void AutoFill2DT<TYPE>::SetHeadRoom(int headroom)
{
	if (headroom < 0) throw eGeneralFail;
	fHeadRoom = headroom;
}

/* accessors */
template <class TYPE>
inline TYPE& AutoFill2DT<TYPE>::operator()(int majordim, int minordim) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
	if (minordim < 0 || minordim >= fCounts[majordim]) throw eOutOfRange;
#endif

	return fArray[majordim*fMaxMinorDim + minordim];
}

template <class TYPE>
inline TYPE* AutoFill2DT<TYPE>::operator()(int majordim) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	return fArray + majordim*fMaxMinorDim;
}

/* set logical size to 0 */
template <class TYPE>
//inline
void AutoFill2DT<TYPE>::Reset(void)
{
	memset(fCounts, 0, sizeof(int)*fMajorDim);
}

template <class TYPE>
inline void AutoFill2DT<TYPE>::Reset(int majordim)
{
#if __option(extended_errorcheck)
	/* checks */
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	fCounts[majordim] = 0;
}

/* assignment operator - dimensions must match exactly */
template <class TYPE>
void AutoFill2DT<TYPE>::Copy(const AutoFill2DT<TYPE>& source)
{
	/* no copies to self */
	if (this != &source)
	{
		/* dimensions must match */
		if (fMajorDim    != source.fMajorDim ||
		    fMaxMinorDim != source.fMaxMinorDim) throw eSizeMismatch;
		
		/* copy counts */
		memcpy(fCounts, source.fCounts, sizeof(int)*fMajorDim);
		
		/* copy data (whole block) */
		memcpy(fArray, source.fArray, sizeof(TYPE)*fMajorDim*fMaxMinorDim);
	}
}

/* adding values to rows */
template <class TYPE>
inline void AutoFill2DT<TYPE>::Append(int majordim, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	int& count = fCounts[majordim];
	
	/* need more memory */
	if (count == fMaxMinorDim) Allocate(fMaxMinorDim + 1, fHeadRoom);
	
	*((*this)(majordim) + count) = value;	
	count++;	
}

template <class TYPE>
void AutoFill2DT<TYPE>::Append(int majordim, const ArrayT<TYPE>& source)
{
#if __option(extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	int length = source.Length();
	if (length > 0)
	{
		int& count  = fCounts[majordim];
		
		/* need more memory */
		if (count + length >= fMaxMinorDim) Allocate(fMaxMinorDim + length, fHeadRoom);
		
		/* copy data */
		memcpy((*this)(majordim) + count, source.Pointer(), sizeof(TYPE)*length);
		count += length;
	}
}

/* add only if not already in the list - returns 1 if */
template <class TYPE>
int AutoFill2DT<TYPE>::AppendUnique(int majordim, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (majordim < 0 || majordim >= fMajorDim) throw eOutOfRange;
#endif

	/* scan logical size for duplicates */
	TYPE* pthis = fArray + majordim*fMaxMinorDim;
	int   count = fCounts[majordim];
	
	for (int i = 0; i < count; i++)
		if (*pthis++ == value)
			return 0;
			
	/* append value on fall through */
	Append(majordim, value);			
	return 1;
}

template <class TYPE>
inline int AutoFill2DT<TYPE>::AppendUnique(int majordim, const ArrayT<TYPE>& source)
{	
	TYPE* psrc = source.Pointer();
	int length = source.Length();
	int count = 0;
	for (int i = 0; i < length; i++)
		count += AppendUnique(majordim, *psrc++);
	
	return count;
}

/***********************************************************************
* Private
***********************************************************************/

template <class TYPE>
void AutoFill2DT<TYPE>::Allocate(int maxminordim, int headroom)
{
	/* initialize count array */
	int copy_in = 1;
	if (!fCounts)
	{
#ifdef __NEW_THROWS__
		try { fCounts = new int[fMajorDim]; }
		catch (bad_alloc) { fCounts = NULL; }
#else
		fCounts = new int[fMajorDim];
#endif

		if (!fCounts)
		{
			cout << "\n AutoFill2DT<TYPE>::Allocate: out of memory"<< endl;
			throw eOutOfMemory;
		}

		memset(fCounts, 0, sizeof(int)*fMajorDim);
		copy_in = 0;
	}

	/* determine new memory size */
	int memsize;
	if (headroom > 0)
		memsize  = (maxminordim*(100 + headroom))/100;
	else
		memsize = maxminordim;
	
	TYPE* newfArray;
#ifdef __NEW_THROWS__
	try { newfArray = new TYPE[fMajorDim*memsize]; }
	catch (bad_alloc) { newfArray = NULL; }
#else
	newfArray = new TYPE[fMajorDim*memsize];
#endif

	if (!newfArray)
	{
		cout << "\n AutoFill2DT<TYPE>::Allocate: out of memory"<< endl;
		throw eOutOfMemory;
	}

	/* keep current data */
	if (copy_in)
	{
		int* pcount = fCounts;
		TYPE* pdata = fArray;
		TYPE* pnew  = newfArray;
	
		for (int i = 0; i < fMajorDim; i++)
		{
			memcpy(pnew, pdata, sizeof(TYPE)*(*pcount++));
			
			pdata += fMaxMinorDim;
			pnew  += memsize;
		}
	}

	delete[] fArray;
	fArray = newfArray;

	fMaxMinorDim = memsize;
}

#endif /* _AUTO_ARRAY2D_T_H_ */
