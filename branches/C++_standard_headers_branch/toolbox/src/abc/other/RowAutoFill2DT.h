/* $Id: RowAutoFill2DT.h,v 1.10 2007-09-05 00:22:47 paklein Exp $ */

#ifndef _ROW_AUTO_ARRAY2D_T_H_
#define _ROW_AUTO_ARRAY2D_T_H_

#include <string.h>
#include <fstream.h>

#include "Environment.h"
#include "ExceptionT.h"

namespace Tahoe {

/** class to allow dynamical resizing of rows. Allocation is on a
 * row by row basis. "Overallocation" is controlled by the headroom
 * parameter. */
template <class TYPE>
class RowAutoFill2DT
{
public:

	/** constructors */
	RowAutoFill2DT(int major_dim, int head_room);
	RowAutoFill2DT(int major_dim, int head_room, int init_row_memory);

	/** destructor */
	~RowAutoFill2DT(void);

	/* dimensions */
	int MajorDim(void) const;
	int HeadRoom(void) const;
	int MinorDim(int major_dim) const;
	int LogicalSize(void) const;
	int MinorDimCount(int dim) const; // returns number of rows with the specified size

	/** re-sizing */
	void SetHeadRoom(int head_room);
	
	/** flush size */
	void SetFlushSize(long flush_size);	
	
	/** \name accessors */
	/*@{*/
	TYPE& operator()(int major_dim, int minor_dim);
	const TYPE& operator()(int major_dim, int minor_dim) const;
	TYPE* operator()(int major_dim);
	const TYPE* operator()(int major_dim) const;
	/*@}*/

	/* set logical size = 0 */
	void Reset(void);    // for all rows
	void Reset(int major_dim); // for the selected row
	
	/** adding value to row */
	void Append(int row, const TYPE& value);

	/** adding source array to row */
	void Append(int row, const ArrayT<TYPE>& source);
	
	/** appending unique value to row. 
	 * \return 1 if value appended, 0 otherwise */
	int AppendUnique(int row, const TYPE& value);

	/** appending unique values from source to row. 
	 * \return number of values appended */
	int AppendUnique(int row, const ArrayT<TYPE>& source);
	
private:
	
	/** set logical size of a row. */
	void SetLogicalSize(int row, int length);

	/** no assigment operator */
	RowAutoFill2DT& operator=(RowAutoFill2DT&);
	
	/** flush memory. write all data to a disk file, free all memory,
	 * reallocate, and read from disk */
	void FlushMemory(void);
	 
private:

	/* dimensions */
	int fHeadRoom; /**< amount of overallocation, as a percentage */
	int fTotalMemorySize; /**< running count of total allocation */
	long fFlushSize; /**< memory allocation disk flush interval in bytes */
	long fFlushCount; /**< allocation count */

	/** logical size of each row */
	ArrayT<int> fLogicalSize;

	/** memory size of each row */
	ArrayT<int> fMemorySize;

	/** row pointers */
	ArrayT<TYPE*> fRowData;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
inline RowAutoFill2DT<TYPE>::RowAutoFill2DT(int major_dim, int head_room):
	fHeadRoom(head_room),
	fTotalMemorySize(0),
	fFlushSize(-1),
	fFlushCount(0),
	fLogicalSize(major_dim),
	fMemorySize(major_dim),
	fRowData(major_dim)
{
	/* check */
	if (fHeadRoom < 0) ExceptionT::GeneralFail("RowAutoFill2DT");
	
	/* initialize sizes */
	fLogicalSize = 0;
	fMemorySize = 0;
	
	/* initialize all pointers */
	fRowData = NULL;
}

template <class TYPE>
RowAutoFill2DT<TYPE>::RowAutoFill2DT(int major_dim, int head_room, int init_row_memory):
	fHeadRoom(head_room),
	fTotalMemorySize(0),
	fFlushSize(-1),
	fFlushCount(0),
	fLogicalSize(major_dim),
	fMemorySize(major_dim),
	fRowData(major_dim)
{
	/* check */
	if (fHeadRoom < 0 || init_row_memory < 0) ExceptionT::GeneralFail("RowAutoFill2DT");

	/* initialize sizes */
	fLogicalSize = 0;
	fMemorySize = 0;
	fRowData = NULL;
	
	/* initialize memory */
	if (init_row_memory > 0)
	{
		/* assume no additional headroom */
		fHeadRoom = 0;
		for (int i = 0; i < major_dim; i++)
			SetLogicalSize(i, init_row_memory);
			
		/* reset logical size */
		fLogicalSize = 0;
	
		/* reset headroom */
		fHeadRoom = head_room;
	}
}

template <class TYPE>
RowAutoFill2DT<TYPE>::~RowAutoFill2DT(void)
{
	/* free pointers */
	int major_dim = fRowData.Length();
	for (int i = 0; i < major_dim; i++)
		delete[] fRowData[i];
}

/* dimensions */
template <class TYPE>
inline int RowAutoFill2DT<TYPE>::MajorDim(void) const { return fRowData.Length(); }

template <class TYPE>
int RowAutoFill2DT<TYPE>::HeadRoom(void) const { return fHeadRoom; }

template <class TYPE>
inline int RowAutoFill2DT<TYPE>::MinorDim(int major_dim) const
{
#if __option(extended_errorcheck)
	/* range check */
	if (major_dim < 0 || major_dim >= MajorDim()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	return fLogicalSize[major_dim];
}

template <class TYPE>
inline int RowAutoFill2DT<TYPE>::LogicalSize(void) const
{
	int sum = 0;
	for (int i = 0; i < fRowData.Length(); i++)
		sum += fLogicalSize[i];
	return sum;
}

template <class TYPE>
inline void RowAutoFill2DT<TYPE>::SetHeadRoom(int head_room)
{
	if (head_room < 0) ExceptionT::GeneralFail("RowAutoFill2DT");
	fHeadRoom = head_room;
}

/* flush size */
template <class TYPE>
inline void RowAutoFill2DT<TYPE>::SetFlushSize(long flush_size)
{
	fFlushSize = flush_size;
	fFlushCount = 0;
}

/* accessors */
template <class TYPE>
inline const TYPE& RowAutoFill2DT<TYPE>::operator()(int major_dim, int minor_dim) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
	if (minor_dim < 0 || minor_dim >= fLogicalSize[major_dim]) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	return *(fRowData[major_dim] + minor_dim);
}

template <class TYPE>
inline TYPE& RowAutoFill2DT<TYPE>::operator()(int major_dim, int minor_dim)
{
#if __option(extended_errorcheck)
	/* checks */
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
	if (minor_dim < 0 || minor_dim >= fLogicalSize[major_dim]) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	return *(fRowData[major_dim] + minor_dim);
}

template <class TYPE>
inline const TYPE* RowAutoFill2DT<TYPE>::operator()(int major_dim) const
{
#if __option(extended_errorcheck)
	/* checks */
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	return fRowData[major_dim];
}

template <class TYPE>
inline TYPE* RowAutoFill2DT<TYPE>::operator()(int major_dim)
{
#if __option(extended_errorcheck)
	/* checks */
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	return fRowData[major_dim];
}

/* set logical size to 0 */
template <class TYPE>
void RowAutoFill2DT<TYPE>::Reset(void)
{
	fLogicalSize = 0;
}

template <class TYPE>
inline void RowAutoFill2DT<TYPE>::Reset(int major_dim)
{
#if __option(extended_errorcheck)
	/* checks */
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	fLogicalSize[major_dim] = 0;
}

/* adding values to rows */
template <class TYPE>
inline void RowAutoFill2DT<TYPE>::Append(int major_dim, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	int size = fLogicalSize[major_dim];
	
	/* dimension row */
	SetLogicalSize(major_dim, size + 1);
	
	/* add value */
	*(fRowData[major_dim] + size) = value;
}

template <class TYPE>
void RowAutoFill2DT<TYPE>::Append(int major_dim, const ArrayT<TYPE>& source)
{
#if __option(extended_errorcheck)
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	int size = fLogicalSize[major_dim];
	
	/* dimension row */
	int add_length = source.Length();
	SetLogicalSize(major_dim, size + add_length);
	
	/* add values */
	TYPE* p = fRowData[major_dim] + size;
	for (int i = 0; i < add_length; i++)
		*p++ = source[i];
}

/* add only if not already in the list - returns 1 if */
template <class TYPE>
int RowAutoFill2DT<TYPE>::AppendUnique(int major_dim, const TYPE& value)
{
#if __option(extended_errorcheck)
	if (major_dim < 0 || major_dim >= fRowData.Length()) ExceptionT::OutOfRange("RowAutoFill2DT");
#endif

	/* scan logical size for duplicates */
	TYPE* pthis = fRowData[major_dim];
	int   count = fLogicalSize[major_dim];
	for (int i = 0; i < count; i++)
		if (*pthis++ == value)
			return 0;
			
	/* append value on fall through */
	Append(major_dim, value);			
	return 1;
}

template <class TYPE>
inline int RowAutoFill2DT<TYPE>::AppendUnique(int major_dim, const ArrayT<TYPE>& source)
{	
	const TYPE* psrc = source.Pointer();
	int length = source.Length();
	int count = 0;
	for (int i = 0; i < length; i++)
		count += AppendUnique(major_dim, *psrc++);
	
	return count;
}

/***********************************************************************
* Private
***********************************************************************/

/* set logical size of a row */
template <class TYPE>
void RowAutoFill2DT<TYPE>::SetLogicalSize(int row, int length)
{
	/* flush */
	if (fFlushSize > 0 && fFlushCount > fFlushSize) FlushMemory();

	/* need more space? */
	if (length > fMemorySize[row])
	{
		/* current size */
		int old_size = fMemorySize[row];
	
		/* new memory size */
		int mem_size = (fHeadRoom > 0) ? (length*(100 + fHeadRoom))/100 : length;
	
		/* allocate */
		TYPE* new_array;
#ifdef __NEW_THROWS__
		try { new_array = new TYPE[mem_size]; }
		catch (bad_alloc) { new_array = NULL; }
#else
		new_array = new TYPE[mem_size];
#endif
		if (!new_array)
			ExceptionT::OutOfMemory("RowAutoFill2DT::SetLogicalSize");
		fMemorySize[row] = mem_size;
		fTotalMemorySize += (mem_size - old_size);
		fFlushCount += mem_size*sizeof(TYPE);

		/* copy old data */
		if (fLogicalSize[row] > 0 && fRowData[row] != NULL) 
			memcpy(new_array, fRowData[row], sizeof(TYPE)*fLogicalSize[row]);

		/* free old memory */
		delete[] fRowData[row];
		
		/* reset pointer */
		fRowData[row] = new_array;
	}
	
	/* set logical size */
	fLogicalSize[row] = length;
}

/* flush memory */
template <class TYPE>
void RowAutoFill2DT<TYPE>::FlushMemory(void)
{
	const char file[] = "RowAutoFill2DT.tmp";
	cout << "\n RowAutoFill2DT<TYPE>::FlushMemory: flushing memory to disk file: \"" 
	     << file << '\"' << endl; 
	ArrayT<TYPE> dump;

	/* dump data */
	ofstream out(file);
	for (int i = 0; i < fRowData.Length(); i++)
	{
		dump.Set(fLogicalSize[i], fRowData[i]);
		dump.WriteBinary(out);
	}
	out.close();

	/* free data memory */
	for (int j = 0; j < fRowData.Length(); j++)
	  {
		delete[] fRowData[j];
		fRowData[j] = NULL;
	  }
	fMemorySize = 0;
	fTotalMemorySize = 0;

	/* re-allocate/read from file */
	int flush_size = fFlushSize;
	fFlushSize = -1;
	ifstream in(file);
	for (int k = 0; k < fRowData.Length(); k++)
	{
	  if (fLogicalSize[k] > 0)
		{
		/* set logical size */
		SetLogicalSize(k, fLogicalSize[k]);
		
		/* read */
		dump.Set(fLogicalSize[k], fRowData[k]);
		dump.ReadBinary(in);
		}
	}
	in.close();

	/* reset */
	fFlushSize = flush_size;
	fFlushCount = 0;

	cout << "\n RowAutoFill2DT<TYPE>::FlushMemory: flushing memory: DONE" << endl; 
}

} // namespace Tahoe 
#endif /* _ROW_AUTO_ARRAY2D_T_H_ */
